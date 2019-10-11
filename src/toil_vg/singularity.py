"""
    Module for running Docker images with Singularity.   
    Derived from https://github.com/BD2KGenomics/toil/blob/master/src/toil/lib/docker.py

    Contains two user-facing functions: singularityCall and singularityCheckOutput
    Example of using singularityCall in a Toil pipeline to index a FASTA file with SAMtools:
        def toil_job(job):
            work_dir = job.fileStore.getLocalTempDir()
            path = job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta')
            parameters = ['faidx', path]
            singularityCall(job, tool='quay.io/ucgc_cgl/samtools:latest', work_dir=work_dir, parameters=parameters)
"""
import base64
import logging
import subprocess
import pipes
import os
import time

_logger = logging.getLogger(__name__)


def singularityCall(job,
               tool,
               parameters=None,
               workDir=None,
               singularityParameters=None,
               outfile=None):
    """
    Throws CalledProcessorError if the Singularity invocation returns a non-zero exit code
    This function blocks until the subprocess call to Singularity returns
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools:latest).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-v`. Destination convention is /data
    :param list[str] singularityParameters: Parameters to pass to Singularity. Default parameters are `--rm`,
            `--log-driver none`, and the mountpoint `-v work_dir:/data` where /data is the destination convention.
             These defaults are removed if singularity_parmaters is passed, so be sure to pass them if they are desired.
    :param file outfile: Pipe output of Singularity call to file handle
    """
    return _singularity(job, tool=tool, parameters=parameters, workDir=workDir, singularityParameters=singularityParameters,
                        outfile=outfile, checkOutput=False)


def singularityCheckOutput(job,
                      tool,
                      parameters=None,
                      workDir=None,
                      singularityParameters=None):
    """
    Returns the stdout from the Singularity invocation (via subprocess.check_output)
    Throws CalledProcessorError if the Singularity invocation returns a non-zero exit code
    This function blocks until the subprocess call to Singularity returns
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools:latest).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-v`. Destination convention is /data
    :param list[str] singularityParameters: Parameters to pass to Singularity. Default parameters are `--rm`,
            `--log-driver none`, and the mountpoint `-v work_dir:/data` where /data is the destination convention.
             These defaults are removed if singularity_parmaters is passed, so be sure to pass them if they are desired.
    :returns: Stdout from the singularity call
    :rtype: str
    """
    return _singularity(job, tool=tool, parameters=parameters, workDir=workDir,
                   singularityParameters=singularityParameters, checkOutput=True)


def _singularity(job,
            tool,
            parameters=None,
            workDir=None,
            singularityParameters=None,
            outfile=None,
            checkOutput=False):
    """
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `--bind`. Destination convention is /data
    :param list[str] singularityrParameters: Parameters to pass to Singularity. Default parameters are the mountpoint
             `--bind work_dir:/data` where /data is the destination convention.
             These defaults are removed if singularity_parmaters is passed, so be sure to pass them if they are desired.
    :param file outfile: Pipe output of Singularity call to file handle
    :param bool checkOutput: When True, this function returns singularity's output.
    """
    if parameters is None:
        parameters = []
    if workDir is None:
        workDir = os.getcwd()

    # Setup the outgoing subprocess call for singularity
    baseSingularityCall = ['singularity', '-q', 'exec']
    if singularityParameters:
        baseSingularityCall += singularityParameters
    else:
        baseSingularityCall += ['-H', '{}:{}'.format(os.path.abspath(workDir), os.environ.get('HOME')), '--pwd', os.environ.get('HOME')]
        
        
    # Problem: Multiple Singularity calls sharing the same cache directory will
    # not work correctly. See https://github.com/sylabs/singularity/issues/3634
    # and https://github.com/sylabs/singularity/issues/4555.
    
    # As a hacky workaround, we use a fresh cache for every Singularity run.
    # TODO: use one per job or per worker.
    singularity_env = os.environ.copy()
    singularity_env['SINGULARITY_CACHEDIR'] = job.fileStore.getLocalTempDir()

    # Make subprocess call

    # If parameters is list of lists, treat each list as separate command and chain with pipes
    if len(parameters) > 0 and type(parameters[0]) is list:
        # When piping, all arguments now get merged into a single string to bash -c.
        # We try to support spaces in paths by wrapping them all in quotes first.
        chain_params = [' '.join(p) for p in [map(pipes.quote, q) for q in parameters]]
        # Use bash's set -eo pipefail to detect and abort on a failure in any command in the chain
        call = baseSingularityCall + [_convertImageSpec(tool), '/bin/bash', '-c',
                                 'set -eo pipefail && {}'.format(' | '.join(chain_params))]
    else:
        call = baseSingularityCall + [_convertImageSpec(tool)] + parameters
    _logger.info("Calling singularity with " + repr(call))

    params = {'env': singularity_env}
    if outfile:
        params['stdout'] = outfile
    if checkOutput:
        callMethod = subprocess.check_output
    else:
        callMethod = subprocess.check_call

    out = callMethod(call, **params)

    # After Singularity exits, it is possible that cleanup of the container's
    # temporary files is still in progress (sicne it also waits for the
    # container to exit). If we return immediately and the Toil job then
    # immediately finishes, we can have a race between Toil's temp
    # cleanup/space tracking code and Singularity's temp cleanup code to delete
    # the same directory tree. Toil doesn't handle this well, and crashes when
    # files it expected to be able to see are missing (at least on some
    # versions). So we introduce a delay here to try and make sure that
    # Singularity wins the race with high probability.
    #
    # See https://github.com/sylabs/singularity/issues/1255
    time.sleep(0.5)

    return out
    
def _convertImageSpec(spec):
    """
    Given an image specifier that may be either a Docker container specifier,
    or a Singularity URL or filename, produce the Singularity URL or filename
    that points to it.
    
    This consists of identifying the Docker container specifiers and prefixing
    them with "docker://".
    """
   
    if spec.startswith('/'):
        # It's a file path we can use.
        # Relative paths won't work because Toil uses unique working
        # directories.
        return spec
   
    if '://' in spec:
        # Already a URL
        return spec
    
    # Try it as a Docker specifier
    return 'docker://' + spec
    
    
