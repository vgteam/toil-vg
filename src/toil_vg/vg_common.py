#!/usr/bin/env python2.7
"""
Shared stuff between different modules in this package.  Some
may eventually move to or be replaced by stuff in toil-lib.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno
from uuid import uuid4

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_lib.programs import docker_call

def add_docker_tool_parse_args(parser):
    """ centralize shared docker options and their defaults """
    parser.add_argument("--no_docker", action="store_true",
                        help="do not use docker for any commands")
    parser.add_argument("--vg_docker", type=str, default='quay.io/glennhickey/vg:latest',
                        help="dockerfile to use for vg")
    parser.add_argument("--bcftools_docker", type=str, default='quay.io/cmarkello/bcftools',
                        help="dockerfile to use for bcftools")
    parser.add_argument("--tabix_docker", type=str, default='quay.io/cmarkello/htslib:latest',
                        help="dockerfile to use for tabix")
    parser.add_argument("--jq_docker", type=str, default='devorbitus/ubuntu-bash-jq-curl',
                        help="dockerfile to use for jq")

def get_docker_tool_map(options):
    """ convenience function to parse the above _docker options into a dictionary """

    dmap = dict()
    if not options.no_docker:
        dmap["vg"] = options.vg_docker
        dmap["bcftools"] = options.bcftools_docker
        dmap["tabix"] = options.tabix_docker
        dmap["bgzip"] = options.tabix_docker
        dmap["jq"] = options.jq_docker

    # to do: could be a good place to do an existence check on these tools

    return dmap
        
class DockerRunner(object):
    """ Helper class to centralize docker calling.  So we can toggle Docker
on and off in just one place.  to do: Should go somewhere more central """
    def __init__(self, docker_tool_map = {}):
        # this maps a command to its full docker name
        # example:  docker_tool_map['vg'] = 'quay.io/ucsc_cgl/vg:latest'
        self.docker_tool_map = docker_tool_map

    def call(self, args, work_dir = '.' , outfile = None, errfile = None,
             check_output = False, inputs=[]):
        """ run a command.  decide to use docker based on whether
        its in the docker_tool_map.  args is either the usual argument list,
        or a list of lists (in the case of a chain of piped commands)  """
        # from here on, we assume our args is a list of lists
        if len(args) == 0 or len(args) > 0 and type(args[0]) is not list:
            args = [args]
        if args[0][0] in self.docker_tool_map:
            return self.call_with_docker(args, work_dir, outfile, errfile, check_output, inputs)
        else:
            return self.call_directly(args, work_dir, outfile, errfile, check_output, inputs)
        
    def call_with_docker(self, args, work_dir, outfile, errfile, check_output, inputs): 
        """ Thin wrapper for docker_call that will use internal lookup to
        figure out the location of the docker file.  Only exposes docker_call
        parameters used so far.  expect args as list of lists.  if (toplevel)
        list has size > 1, then piping interface used """

        RealTimeLogger.get().info("Docker Run: {}".format(" | ".join(" ".join(x) for x in args)))

        if len(args) == 1:
            # just one command, use regular docker_call interface
            # where parameters is an argument list not including command
            tool = self.docker_tool_map[args[0][0]]
            tools = None
            if args[0][0] == "vg":
                # todo:  this is a hack because the vg docker currently *does not* expect
                # command lines passed in to begin with vg.  Seems inconsistent with
                # all other containers (ex bcftools expects bcftools as first arg)
                # it's very hard to work around programmatically when supporting
                # docker/non-docker consistency to have to parse different command lines
                # differently. 
                parameters = args[0][1:]
            else:
                parameters = args[0]
        else:
            # there's a pipe.  we use the different piping interface that
            # takes in paramters as a list of single-string commands
            # that include arguments
            tool = None
            tools = self.docker_tool_map[args[0][0]]
            parameters = [" ".join(x) for x in args]

        return docker_call(tool=tool, tools=tools, parameters=parameters,
                           work_dir=work_dir, outfile = outfile,
                           errfile = errfile,
                           check_output = check_output,
                           inputs=inputs)

    def call_directly(self, args, work_dir, outfile, errfile, check_output, inputs):
        """ Just run the command without docker """

        RealTimeLogger.get().info("Run: {}".format(" | ".join(" ".join(x) for x in args)))

        # this is all that docker_call does with the inputs parameter:
        for filename in inputs:
            assert(os.path.isfile(os.path.join(work_dir, filename)))

        procs = []
        for i in range(len(args)):
            stdin = procs[i-1].stdout if i > 0 else None
            if i == len(args) - 1 and outfile is not None:
                stdout = outfile
            else:
                stdout = subprocess.PIPE
            procs.append(subprocess.Popen(args[i], stdout=stdout, stderr=errfile,
                                          stdin=stdin, cwd=work_dir))
            
        for p in procs[:-1]:
            p.stdout.close()

        output, errors = procs[-1].communicate()
        for i, proc in enumerate(procs):
            sts = proc.wait()
            if sts != 0:            
                raise Exception("Command {} returned with non-zero exit status {}".format(
                    " ".join(args[i]), sts))

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

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    
    From http://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch
