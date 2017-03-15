"""
    Module for calling Singularity. Assumes `singularity` is on the PATH.
    Contains two user-facing functions: singularityCall and singularityCheckOutput
    Example of using singularityCall in a Toil pipeline to index a FASTA file with SAMtools:
        def toil_job(job):
            work_dir = job.fileStore.getLocalTempDir()
            path = job.fileStore.readGlobalFile(ref_id, os.path.join(work_dir, 'ref.fasta')
            parameters = ['faidx', path]
            dockerCall(job, tool='quay.io/ucgc_cgl/samtools:latest', work_dir=work_dir, parameters=parameters)
"""
import base64
import logging
import subprocess
import pipes
import os
from bd2k.util.exceptions import require

_logger = logging.getLogger(__name__)


def singularityCall(job,
               tool,
               parameters=None,
               workDir=None,
               singularityParameters=None,
               outfile=None):
    """
    Throws CalledProcessorError if the Docker invocation returns a non-zero exit code
    This function blocks until the subprocess call to Docker returns
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Singularity image to be used (e.g. quay.io/ucsc_cgl/samtools:latest).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-v`. Destination convention is /data
    :param list[str] dockerParameters: Parameters to pass to Docker. Default parameters are `--rm`,
            `--log-driver none`, and the mountpoint `-v work_dir:/data` where /data is the destination convention.
             These defaults are removed if docker_parmaters is passed, so be sure to pass them if they are desired.
    :param file outfile: Pipe output of Docker call to file handle
    """
    _singularity(job, tool=tool, parameters=parameters, workDir=workDir, singularityParameters=singularityParameters,
            outfile=outfile, checkOutput=False, defer=defer)


def singularityCheckOutput(job,
                      tool,
                      parameters=None,
                      workDir=None,
                      singularityParameters=None):
    """
    Returns the stdout from the Docker invocation (via subprocess.check_output)
    Throws CalledProcessorError if the Docker invocation returns a non-zero exit code
    This function blocks until the subprocess call to Docker returns
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools:latest).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-v`. Destination convention is /data
    :param list[str] dockerParameters: Parameters to pass to Docker. Default parameters are `--rm`,
            `--log-driver none`, and the mountpoint `-v work_dir:/data` where /data is the destination convention.
             These defaults are removed if docker_parmaters is passed, so be sure to pass them if they are desired.
    :returns: Stdout from the docker call
    :rtype: str
    """
    return _singularity(job, tool=tool, parameters=parameters, workDir=workDir,
                   singularityParameters=singularityParameters, checkOutput=True, defer=defer)


def _singularity(job,
            tool,
            parameters=None,
            workDir=None,
            singularityParameters=None,
            outfile=None,
            checkOutput=False,
            defer=None):
    """
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Singularity image to be used (e.g. quay.io/ucsc_cgl/samtools).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `--bind`. Destination convention is /data
    :param list[str] dockerParameters: Parameters to pass to Singularity. Default parameters are the mountpoint
             `--bind work_dir:/data` where /data is the destination convention.
             These defaults are removed if singularity_parmaters is passed, so be sure to pass them if they are desired.
    :param file outfile: Pipe output of Singularity call to file handle
    :param bool checkOutput: When True, this function returns singularity's output.
    """
    if parameters is None:
        parameters = []
    if workDir is None:
        workDir = os.getcwd()

    # Setup the outgoing subprocess call for docker
    baseSingularityCall = ['singularity', 'exec']
    if singularityParameters:
        baseSingularityCall += singularityParameters
    else:
        baseSingularityCall += ['--bind',
                           os.path.abspath(workDir) + ':/data']

    # Make subprocess call

    # If parameters is list of lists, treat each list as separate command and chain with pipes
    if len(parameters) > 0 and type(parameters[0]) is list:
        # When piping, all arguments now get merged into a single string to bash.
        # We try to support spaces in paths by wrapping them all in quotes first.
        chain_params = [' '.join(p) for p in parameters]
        call = baseSingularityCall + [tool, ' | '.join(chain_params)]
    else:
        call = baseSingularityCall + [tool] + parameters
    
    call = " ".join(call)
    _logger.info("Calling singularity with " + repr(call))
    
    if outfile:
        subprocess.check_call(call, stdout=outfile, shell=True)
    else:
        if checkOutput:
            return subprocess.check_output(call, shell=True)
        else:
            subprocess.check_call(call, shell=True)

