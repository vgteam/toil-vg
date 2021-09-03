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
import hashlib
import logging
import subprocess
import pipes
import os
import pathlib
import shutil
import sys
import tempfile
import time

logger = logging.getLogger(__name__)

def is_containerized():
    """
    Return True if we think we are already running in a Docker/Kubernetes
    container (where Singularity is unlikely to work without user-mode
    namespaces), and False otherwsie.
    """

    if not os.path.exists('/proc/self/cgroup'):
        # Not on container-having Linux
        return False

    with open('/proc/self/cgroup') as fh:
        for line in fh:
            line = line.lower()
            if 'docker' in line or 'kube' in line:
                # If any of the cgroups smells Docker or Kube-y, assume we are
                # in a container.
                return True
    return False


def singularityCall(job,
               tool,
               parameters=None,
               workDir=None,
               singularityParameters=None,
               outfile=None,
               mount_list=None):
    """
    Throws CalledProcessorError if the Singularity invocation returns a non-zero exit code
    This function blocks until the subprocess call to Singularity returns
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools:latest).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-B`.
           Destination convention is /mnt, which almost certainly exists in the
           container.
    :param list[str] singularityParameters: Parameters to pass to Singularity.
           Overrides defaults which mount the workDir and configure user mode and
           writability.
    :param file outfile: Pipe output of Singularity call to file handle
    :param list[str] mount_list: List of directories from which to mount into the
           container via `-B`. Destination convention is /mnt, which almost certainly
           exists in the container.
    """
    return _singularity(job, tool=tool, parameters=parameters, workDir=workDir, singularityParameters=singularityParameters,
                        outfile=outfile, checkOutput=False, mount_list=mount_list)


def singularityCheckOutput(job,
                      tool,
                      parameters=None,
                      workDir=None,
                      singularityParameters=None,
                      mount_list=None):
    """
    Returns the stdout from the Singularity invocation (via subprocess.check_output)
    Throws CalledProcessorError if the Singularity invocation returns a non-zero exit code
    This function blocks until the subprocess call to Singularity returns
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools:latest).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-B`.
           Destination convention is /mnt, which almost certainly exists in the
           container.
    :param list[str] singularityParameters: Parameters to pass to Singularity.
           Overrides defaults which mount the workDir and configure user mode and
           writability.
    :param list[str] mount_list: List of directories from which to mount into the
           container via `-B`. Destination convention is /mnt, which almost certainly
           exists in the container.
    :returns: Stdout from the singularity call
    :rtype: str
    """
    return _singularity(job, tool=tool, parameters=parameters, workDir=workDir,
                   singularityParameters=singularityParameters, checkOutput=True, mount_list=mount_list)


def _singularity(job,
            tool,
            parameters=None,
            workDir=None,
            singularityParameters=None,
            outfile=None,
            checkOutput=False,
            mount_list=None):
    """
    :param toil.Job.job job: The Job instance for the calling function.
    :param str tool: Name of the Docker image to be used (e.g. quay.io/ucsc_cgl/samtools).
    :param list[str] parameters: Command line arguments to be passed to the tool.
           If list of lists: list[list[str]], then treat as successive commands chained with pipe.
    :param str workDir: Directory to mount into the container via `-B`.
           Destination convention is /mnt, which almost certainly exists in the
           container.
    :param list[str] singularityParameters: Parameters to pass to Singularity.
           Overrides defaults which mount the workDir and configure user mode and
           writability.
    :param file outfile: Pipe output of Singularity call to file handle
    :param bool checkOutput: When True, this function returns singularity's output.
    :param list[str] mount_list: List of directories from which to mount into the
           container via `-B`. Destination convention is /mnt, which almost certainly
           exists in the container.
    """
    if parameters is None:
        parameters = []
    if workDir is None:
        workDir = os.getcwd()

    # Setup the outgoing subprocess call for singularity.
    baseSingularityCall = ['singularity', 'exec']
    if singularityParameters:
        baseSingularityCall += singularityParameters
    else:
        # Make the container writable. Writing to the container is still not
        # advised, but some tools/environments break if it is not enabled.
        baseSingularityCall.append('-w')

        if is_containerized():
            # We are already in a container. We need to run in user mode,
            # because the container may not be privileged. If we don't use user
            # mode, Singularity tries to do some confining that it can't do in
            # an un-privileged container, and fails.
            baseSingularityCall.append('-u')
            
        if not str(pathlib.Path.home()).startswith('/home'):
            # Newer versions of Singularity will fail if they can't mount the
            # home directory, which they can't do if the directory the home
            # directory is in doesn't exist in the container. If it isn't just
            # under /home, assume it might not be in the container and tell
            # Singularity not to try to mount it. See
            # https://github.com/hpcng/singularity/issues/4995
            baseSingularityCall.append('--no-home')

        # Mount workdir as /mnt and work in there.
        # Hope the image actually has a /mnt available.
        # Otherwise this silently doesn't mount.
        # But with -u (user namespaces) we have no luck pointing in-container
        # home at anything other than our real home (like something under /var
        # where Toil puts things).
        # Note that we target Singularity 3+.
        mount_string =  '{}:{}'.format(os.path.abspath(workDir), '/mnt')
        if mount_list:
            for mount_dir in mount_list:
                mount_string = '{},{}:{}'.format(mount_string, os.path.abspath(mount_dir), '/mnt/{}'.format(os.path.basename(mount_dir)))
        
        baseSingularityCall += ['-B', mount_string, '--pwd', '/mnt']
        
    # Problem: Multiple Singularity downloads sharing the same cache directory will
    # not work correctly. See https://github.com/sylabs/singularity/issues/3634
    # and https://github.com/sylabs/singularity/issues/4555.
    
    # As a workaround, we have out own cache which we manage ourselves.
    cache_dir = os.path.join(os.environ.get('SINGULARITY_CACHEDIR',  os.path.join(os.environ.get('HOME'), '.singularity')), 'toil')
    os.makedirs(cache_dir, exist_ok=True)
    
    # What Singularity url/spec do we want?
    source_image = _convertImageSpec(tool) 
    
    # What name in the cache dir do we want?
    # We cache everything as sandbox directories and not .sif files because, as
    # laid out in https://github.com/sylabs/singularity/issues/4617, there
    # isn't a way to run from a .sif file and have write permissions on system
    # directories in the container, because the .sif build process makes
    # everything owned by root inside the image. Since some toil-vg containers
    # (like the R one) want to touch system files (to install R packages at
    # runtime), we do it this way to act more like Docker.
    #
    # Also, only sandbox directories work with user namespaces, and only user
    # namespaces work inside unprivileged Docker containers like the Toil
    # appliance.
    sandbox_dirname = os.path.join(cache_dir, '{}.sandbox'.format(hashlib.sha256(source_image.encode()).hexdigest()))
    
    if not os.path.exists(sandbox_dirname):
        # We atomically drop the sandbox at that name when we get it
        
        # Make a temp directory to be the sandbox
        temp_sandbox_dirname = tempfile.mkdtemp(dir=cache_dir)

        # Download with a fresh cache to a sandbox
        download_env = os.environ.copy()
        download_env['SINGULARITY_CACHEDIR'] = job.fileStore.getLocalTempDir()
        subprocess.check_call(['singularity', 'build', '-s', '-F', temp_sandbox_dirname, source_image], env=download_env)
        
        # Clean up the Singularity cache since it is single use
        shutil.rmtree(download_env['SINGULARITY_CACHEDIR'])
        
        try:
            # This may happen repeatedly but it is atomic
            os.rename(temp_sandbox_dirname, sandbox_dirname)
        except (FileExistsError, OSError) as e:
            # Can't rename a directory over another
            # Make sure someone else has made the directory
            assert os.path.exists(sandbox_dirname)
            # Remove our redundant copy
            shutil.rmtree(temp_sandbox_dirname)
            
        # TODO: we could save some downloading by having one process download
        # and the others wait, but then we would need a real fnctl locking
        # system here.

    # Make subprocess call for singularity run
    
    # Set the TMPDIR environment variable to relative path '.' when running an already built
    # container. This is to get around issues with differences between TMPDIR path accessibility
    # when running Singularity build and singularity exec
    download_env = os.environ.copy()
    if not 'rocker/tidyverse' in tool: 
        download_env['TMPDIR'] = '.'
    
    # If parameters is list of lists, treat each list as separate command and chain with pipes
    if len(parameters) > 0 and type(parameters[0]) is list:
        # When piping, all arguments now get merged into a single string to bash -c.
        # We try to support spaces in paths by wrapping them all in quotes first.
        chain_params = [' '.join(p) for p in [list(map(pipes.quote, q)) for q in parameters]]
        # Use bash's set -eo pipefail to detect and abort on a failure in any command in the chain
        call = baseSingularityCall + [sandbox_dirname, '/bin/bash', '-c',
                                 'set -eo pipefail && {}'.format(' | '.join(chain_params))]
    else:
        call = baseSingularityCall + [sandbox_dirname] + parameters
    logger.info("Calling singularity with " + repr(call))

    params = {}
    params['env'] = download_env
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
    
    
