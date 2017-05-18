#!/usr/bin/env python2.7
"""
Shared stuff between different modules in this package.  Some
may eventually move to or be replaced by stuff in toil-lib.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno
from uuid import uuid4
import pkg_resources, tempfile, datetime
import logging

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.docker import dockerCall, dockerCheckOutput, _fixPermissions
from toil_vg.singularity import singularityCall, singularityCheckOutput
from toil_vg.iostore import IOStore

logger = logging.getLogger(__name__)

def add_container_tool_parse_args(parser):
    """ centralize shared container options and their defaults """

    parser.add_argument("--vg_docker", type=str,
                        help="Docker image to use for vg")
    parser.add_argument("--container", default=None, choices=['Docker', 'Singularity', 'None'],
                       help="Container type used for running commands. Use None to "
                       " run locally on command line")    

def add_common_vg_parse_args(parser):
    """ centralize some shared io functions and their defaults """
    parser.add_argument('--config', default=None, type=str,
                        help='Config file.  Use toil-vg generate-config to see defaults/create new file')
    
    parser.add_argument("--force_outstore", action="store_true",
                        help="use output store instead of toil for all intermediate files (use only for debugging)")
                        
    
def get_container_tool_map(options):
    """ convenience function to parse the above _container options into a dictionary """

    dmap = dict()
    smap = dict()
    # first check for docker use, if docker use is not desired then check
    # for singularity use
    if options.container == 'Docker':
        dmap["vg"] = options.vg_docker
        dmap["bcftools"] = options.bcftools_docker
        dmap["tabix"] = options.tabix_docker
        dmap["bgzip"] = options.tabix_docker
        dmap["jq"] = options.jq_docker
        dmap["rtg"] = options.rtg_docker
        dmap["pigz"] = options.pigz_docker
        dmap["samtools"] = options.samtools_docker
        dmap["bwa"] = options.bwa_docker
    elif options.container == 'Singularity':
        smap["vg"] = options.vg_singularity
        smap["bcftools"] = options.bcftools_singularity
        smap["tabix"] = options.tabix_singularity
        smap["bgzip"] = options.tabix_singularity
        smap["jq"] = options.jq_singularity
        smap["rtg"] = options.rtg_singularity
        smap["pigz"] = options.pigz_singularity
        smap["samtools"] = options.samtools_singularity
        smap["bwa"] = options.bwa_singularity
        
    # to do: could be a good place to do an existence check on these tools

    return dmap, smap
        
class ContainerRunner(object):
    """ Helper class to centralize container calling.  So we can toggle both
Docker and Singularity on and off in just one place.
to do: Should go somewhere more central """
    def __init__(self, container_tool_map = {}):
        # this maps a command to its full docker name or singularity name
        #   the first index is a dictionary containing docker tool names
        #   the second index is a dictionary containing singularity tool names
        # example:  docker_tool_map['vg'] = 'quay.io/ucsc_cgl/vg:latest'
        # example:  singularity_tool_map['vg'] = 'shub://cmarkello/vg'
        self.docker_tool_map = container_tool_map[0]
        self.singularity_tool_map = container_tool_map[1]

    def has_tool(self, tool):
        # return true if we have an image for this tool
        return tool in self.docker_tool_map or tool in self.singularity_tool_map

    def call(self, job, args, work_dir = '.' , outfile = None, errfile = None,
             check_output = False, tool_name=None):
        """ run a command.  decide to use docker based on whether
        its in the docker_tool_map.  args is either the usual argument list,
        or a list of lists (in the case of a chain of piped commands)  """
        # from here on, we assume our args is a list of lists
        if len(args) == 0 or len(args) > 0 and type(args[0]) is not list:
            args = [args]
        # convert everything to string
        for i in range(len(args)):
            args[i] = [str(x) for x in args[i]]
        name = tool_name if tool_name is not None else args[0][0]

        if name in self.docker_tool_map and self.docker_tool_map[name] and\
           self.docker_tool_map[name].lower() != 'none':
            return self.call_with_docker(job, args, work_dir, outfile, errfile, check_output, tool_name)
        elif name in self.singularity_tool_map:
            return self.call_with_singularity(job, args, work_dir, outfile, errfile, check_output, tool_name)
        else:
            return self.call_directly(args, work_dir, outfile, errfile, check_output)
        
    def call_with_docker(self, job, args, work_dir, outfile, errfile, check_output, tool_name): 
        """ Thin wrapper for docker_call that will use internal lookup to
        figure out the location of the docker file.  Only exposes docker_call
        parameters used so far.  expect args as list of lists.  if (toplevel)
        list has size > 1, then piping interface used """

        RealtimeLogger.info("Docker Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # we use the first argument to look up the tool in the docker map
        # but allow overriding of this with the tool_name parameter
        name = tool_name if tool_name is not None else args[0][0]
        tool = self.docker_tool_map[name]

        # default parameters from toil's docker.py
        docker_parameters = ['--rm', '--log-driver', 'none']

        if len(args) == 1:
            # split off first argument as entrypoint (so we can be oblivious as to whether
            # that happens by default)
            parameters = [] if len(args[0]) == 1 else args[0][1:]
            docker_parameters += ['--entrypoint', args[0][0]]
        else:
            # can leave as is for piped interface which takes list of args lists
            # and doesn't worry about entrypoints since everything goes through bash -c
            # todo: check we have a bash entrypoint!
            parameters = args
        
        # vg uses TMPDIR for temporary files
        # this is particularly important for gcsa, which makes massive files.
        # we will default to keeping these in our working directory
        docker_parameters += ['--env', 'TMPDIR=.']

        # set our working directory map
        if work_dir is not None:
            docker_parameters += ['-v', '{}:/data'.format(os.path.abspath(work_dir)),
                                  '-w', '/data']

        if check_output is True:
            ret = dockerCheckOutput(job, tool, parameters=parameters,
                                    dockerParameters=docker_parameters, workDir=work_dir)
        else:
            ret = dockerCall(job, tool, parameters=parameters, dockerParameters=docker_parameters,
                             workDir=work_dir, outfile = outfile)
        
        # This isn't running through reliably by itself.  Will assume it's
        # because I took docker.py out of toil, and leave here until we revert to
        # toil's docker
        #
        # Note: It's the tabix docker call in merge_vcf_chunks that's problematic
        #       It complains that the container can't be found, so fixPermissions
        #       doesn't get run afterward.  
        #
        _fixPermissions(tool, work_dir)

        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Successfully docker ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))

        return ret
    
    def call_with_singularity(self, job, args, work_dir, outfile, errfile, check_output, tool_name): 
        """ Thin wrapper for singularity_call that will use internal lookup to
        figure out the location of the singularity file.  Only exposes singularity_call
        parameters used so far.  expect args as list of lists.  if (toplevel)
        list has size > 1, then piping interface used """

        RealtimeLogger.info("Singularity Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # we use the first argument to look up the tool in the singularity map
        # but allow overriding of this with the tool_name parameter
        name = tool_name if tool_name is not None else args[0][0]
        smap_val = self.singularity_tool_map[name]
        tool, entrypoint = smap_val[0], smap_val[1]
        
        if len(args) == 1:
            # one command: we check entry point (stripping first arg if necessary)
            # and pass in single args list
            if entrypoint:
                parameters = args[0][1:]
            else:
                parameters = args[0]
        else:
            # can leave as is for piped interface which takes list of args lists
            # and doesn't worry about entrypoints since everything goes through bash
            parameters = args

        singularity_parameters = None
        # vg uses TMPDIR for temporary files
        # this is particularly important for gcsa, which makes massive files.
        # we will default to keeping these in our working directory
        if work_dir is not None:
            singularity_parameters = ['-H', '{}:/data'.format(os.path.abspath(work_dir)), '--bind', '{}:/data'.format(os.path.abspath(work_dir))]

        if check_output is True:
            ret = singularityCheckOutput(job, tool, parameters=parameters,
                                    singularityParameters=singularity_parameters, workDir=work_dir)
        else:
            ret = singularityCall(job, tool, parameters=parameters, singularityParameters=singularity_parameters,
                             workDir=work_dir, outfile = outfile)
        
        # This isn't running through reliably by itself.  Will assume it's
        # because I took docker.py out of toil, and leave here until we revert to
        # toil's docker
        #
        # Note: It's the tabix docker call in merge_vcf_chunks that's problematic
        #       It complains that the container can't be found, so fixPermissions
        #       doesn't get run afterward.  
        #
        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Successfully singularity ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))

        return ret

    def call_directly(self, args, work_dir, outfile, errfile, check_output):
        """ Just run the command without docker """

        RealtimeLogger.info("Run: {}".format(" | ".join(" ".join(x) for x in args)))
        start_time = timeit.default_timer()

        # vg uses TMPDIR for temporary files
        # this is particularly important for gcsa, which makes massive files.
        # we will default to keeping these in our working directory
        my_env = os.environ.copy()
        my_env['TMPDIR'] = '.'

        procs = []
        for i in range(len(args)):
            stdin = procs[i-1].stdout if i > 0 else None
            if i == len(args) - 1 and outfile is not None:
                stdout = outfile
            else:
                stdout = subprocess.PIPE

            procs.append(subprocess.Popen(args[i], stdout=stdout, stderr=errfile,
                                          stdin=stdin, cwd=work_dir, env=my_env))
            
        for p in procs[:-1]:
            p.stdout.close()

        output, errors = procs[-1].communicate()
        for i, proc in enumerate(procs):
            sts = proc.wait()
            if sts != 0:            
                raise Exception("Command {} returned with non-zero exit status {}".format(
                    " ".join(args[i]), sts))

        end_time = timeit.default_timer()
        run_time = end_time - start_time
        RealtimeLogger.info("Successfully ran {} in {} seconds.".format(
            " | ".join(" ".join(x) for x in args), run_time))            

        if check_output:
            return output

def get_files_by_file_size(dirname, reverse=False):
    """ Return list of file paths in directory sorted by file size """

    # Get list of files
    filepaths = []
    for basename in os.listdir(dirname):
        filename = os.path.join(dirname, basename)
        if os.path.isfile(filename):
            filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i in xrange(len(filepaths)):
        filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

    return filepaths

def clean_toil_path(path):
    """ Try to make input path into something toil friendly """
    # local path
    if ':' not in path:
        return 'file://' + os.path.abspath(path)
    else:
        return path

def init_out_store(options, command):
    """
    Write a little bit of logging to the output store.
    
    Rely on IOStore to create the store if it doesn't exist
    as well as to check its a valid location. 

    Do this at very beginning to avoid finding an outstore issue
    after hours spent computing
     
    """
    f = tempfile.NamedTemporaryFile(delete=True)
    now = datetime.datetime.now()
    f.write('{}\ntoil-vg {} version {}\nOptions:'.format(now, command,
                    pkg_resources.get_distribution('toil-vg').version))
    for key,val in options.__dict__.items():
        f.write('{}: {}\n'.format(key, val))
    f.flush()
    IOStore.get(options.out_store).write_output_file(f.name, 'toil-vg-{}.txt'.format(command))
    f.close()

def import_to_store(toil, options, path, use_out_store = None,
                    out_store_key = None):
    """
    Imports a path into the File or IO store

    Abstract all store writing here so we can switch to the out_store
    when we want to checkpoint an intermeidate file for output
    or just have all intermediate files in the outstore for debugging.
    
    Returns the id in job's file store if use_out_store is True
    otherwise a key (up to caller to make sure its unique)
    in the out_store

    By default options.force_outstore is used to toggle between file and 
    i/o store.  This will be over-ridden by the use_out_store parameter 
    if the latter is not None
    """
    logger.info("Importing {}".format(path))

    if use_out_store is True or (use_out_store is None and options.force_outstore is True):
        return  write_to_store(None, options, path, use_out_store, out_store_key)
    else:
        return toil.importFile(clean_toil_path(path))
    
def write_to_store(job, options, path, use_out_store = None,
                   out_store_key = None):
    """
    Writes path into the File or IO store (from options.out_store)

    Abstract all store writing here so we can switch to the out_store
    when we want to checkpoint an intermeidate file for output
    or just have all intermediate files in the outstore for debugging.
    
    Returns the id in job's file store if use_out_store is True
    otherwise a key (up to caller to make sure its unique)
    in the out_store

    By default options.force_outstore is used to toggle between file and 
    i/o store.  This will be over-ridden by the use_out_store parameter 
    if the latter is not None
    """
    if use_out_store is True or (use_out_store is None and options.force_outstore is True):
        out_store = IOStore.get(options.out_store)
        key = os.path.basename(path) if out_store_key is None else out_store_key
        out_store.write_output_file(path, key)
        return key
    else:
        return job.fileStore.writeGlobalFile(path)

def read_from_store(job, options, id_or_key, path = None, use_out_store = None):
    """
    Reads id (or key) from the File store (or IO store from options.out_store) into path

    Abstract all store reading here so we can switch to the out_store
    when we want to checkpoint an intermeidate file for output
    or just have all intermediate files in the outstore for debugging.

    By default options.force_outstore is used to toggle between file and 
    i/o store.  This will be over-ridden by the use_out_store parameter 
    if the latter is not None
    """
    if use_out_store is True or (use_out_store is None and options.force_outstore is True):
        # we can add this interface if we really want by coming up
        # with unique name here, i guess
        assert path is not None
        out_store = IOStore.get(options.out_store)
        return out_store.read_input_file(id_or_key, path)
    else:
        return job.fileStore.readGlobalFile(id_or_key, path)

def write_dir_to_store(job, options, path, use_out_store = None):
    """
    Need directory interface for rocksdb indexes.  Want to avoid overhead
    of tar-ing up as they may be big.  Write individual files instead, and 
    keep track of the names as well as ids (returns list of name/id pairs)
    """
    out_pairs = []
    file_list = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]
    for f in file_list:
        f_id = write_to_store(job, options, f, use_out_store = use_out_store,
                              out_store_key = path.replace('/', '_'))
        out_pairs.append(os.path.basename(f), f_id)
    return out_pairs

def read_dir_from_store(job, options, name_id_pairs, path = None, use_out_store = None):
    """
    Need directory interface for rocksdb indexes.  Want to avoid overhead
    of tar-ing up as they may be big.  Takes as input list of filename/id pairs
    and reads these into the local directory given
    """
    if not os.path.isdir(path):
        os.mkdir(path)

    for name, key in name_id_pairs:
        read_from_store(job, options, key, os.path.join(path, name),
                        use_out_store = use_out_store)
    
def require(expression, message):
    if not expression:
        raise Exception('\n\n' + message + '\n\n')

def parse_id_ranges(job, options, id_ranges_file_id):
    """Returns list of triples chrom, start, end
    """
    work_dir = job.fileStore.getLocalTempDir()
    id_range_file = os.path.join(work_dir, 'id_ranges.tsv')
    read_from_store(job, options, id_ranges_file_id, id_range_file)
    return parse_id_ranges_file(id_range_file)

def parse_id_ranges_file(id_ranges_filename):
    """Returns list of triples chrom, start, end
    """
    id_ranges = []
    with open(id_ranges_filename) as f:
        for line in f:
            toks = line.split()
            if len(toks) == 3:
                id_ranges.append((toks[0], int(toks[1]), int(toks[2])))
    return id_ranges
                                 
                
