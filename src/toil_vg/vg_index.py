#!/usr/bin/env python2.7
"""
vg_index.py: index a graph so it can be mapped to

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import getpass
import pdb

from math import ceil
from subprocess import Popen, PIPE
from Bio import SeqIO

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_lib.programs import docker_call
from toil_vg.vg_common import *

def parse_args():
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """

    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(prog='vg_evaluation_pipeline', description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("vg_graph", type=str,
        help="Input vg graph file path")
    parser.add_argument("out_dir", type=str,
        help="directory where all output will be written")
    parser.add_argument("input_store",
        help="sample input IOStore where input files will be temporarily uploaded")
    parser.add_argument("out_store",
        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")
    parser.add_argument("--path_name", nargs='+', type=str,
        help="Name of reference path in the graph (eg. ref or 17)")    

    # Add common docker options
    add_docker_tool_parse_args(parser)

    # Add indexing options
    index_parse_args(parser)

    options = parser.parse_args()

    return parser.parse_args()


def index_parse_args(parser):
    """ centralize indexing parameters here """
    
    parser.add_argument("--edge_max", type=int, default=5,
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int, default=10,
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--reindex", default=False, action="store_true",
        help="don't re-use existing indexed graphs")
    parser.add_argument("--index_mode", choices=["rocksdb", "gcsa-kmer",
        "gcsa-mem"], default="gcsa-mem",
        help="type of vg index to use for mapping")
    parser.add_argument("--include_pruned", action="store_true",
        help="use the pruned graph in the index")
    parser.add_argument("--include_primary", action="store_true",
        help="use the primary path in the index")
    parser.add_argument("--index_cores", type=int, default=3,
        help="number of threads during the indexing step")

    
def build_indexes(job, options):
    """ Common indexing logic that's shared between standalone indexing and 
    vg_evaluation_pipeline.  This is run from a toil job but not actually 
    added as a child or follow-on in either case.  The point of pulling this
    out is that run_indexing() in its current form has a child that isn't used
    in standalone.  So we re-use this function within different toil hierarchies
    todo: this shouldn't be that hard to clean up """
    
    RealTimeLogger.get().info("Starting indexing...")
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    graph_file = os.path.basename(options.vg_graph)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download local input files from the remote storage container
    graph_dir = work_dir
    robust_makedirs(graph_dir)
    
    graph_filename = "{}/graph.vg".format(graph_dir)
    graph_file_remote_path = graph_file
    input_store.read_input_file(graph_file_remote_path, graph_filename)
    
    
    # Now run the indexer.
    RealTimeLogger.get().info("Indexing {}".format(options.vg_graph))
            
    if options.index_mode == "rocksdb":
        # Make the RocksDB index
        command = ['index', '-s', '-k', str(options.kmer_size), '-e', str(options.edge_max), '-t', str(job.cores), os.path.basename(graph_filename), os.path.basename('{}/{}.index'.format(graph_dir, graph_file))]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest')

    elif (options.index_mode == "gcsa-kmer" or
        options.index_mode == "gcsa-mem"):
        # We want a GCSA2/xg index. We have to prune the graph ourselves.
        # See <https://github.com/vgteam/vg/issues/286>.

        # What will we use as our temp combined graph file (containing only
        # the bits of the graph we want to index, used for deduplication)?
        to_index_filename = "{}/to_index.vg".format(
            work_dir)

        # Where will we save the kmers?
        kmers_filename = "{}/index.graph".format(
            work_dir)

        
        with open(to_index_filename, "w") as to_index_file:

            if options.include_pruned:

                RealTimeLogger.get().info("Pruning {} to {}".format(
                    graph_filename, to_index_filename))

                # Prune out hard bits of the graph
                # and complex regions
                # and short disconnected chunks
                command = ['vg mod -p -l {} -t {} -e {} {}'.format(str(options.kmer_size), str(job.cores), str(options.edge_max), os.path.basename(graph_filename)),
                            'vg mod -S -l {} -t {} -'.format(str(options.kmer_size * 2), str(job.cores))]
                docker_call(work_dir=work_dir, parameters=command,
                            tools='quay.io/ucsc_cgl/vg:latest',
                            outfile=to_index_file)

            if options.include_primary:

                # Then append in the primary path. Since we don't knoiw what
                # "it's called, we retain "ref" and all the 19", "6", etc paths
                # "from 1KG.

                RealTimeLogger.get().info(
                    "Adding primary path to {}".format(to_index_filename))
            
                RealTimeLogger.get().info(
                    "to_index_file: {}".format(to_index_file))
                # See
                # https://github.com/vgteam/vg/issues/318#issuecomment-215102199

                # Generate all the paths names we might have for primary paths.
                # It should be "ref" but some graphs don't listen
                RealTimeLogger.get().info("OPTIONS.PATH_NAME: {}".format(options.path_name[0]))
                RealTimeLogger.get().info("CHR IN OPTIONS.PATH_NAME: {}".format('chr' in options.path_name[0]))
                if 'chr' in options.path_name[0]:
                    ref_names = (["ref", "chrx", "chrX", "chry", "chrY", "chrm", "chrM"] +
                        ['chr'+str(x) for x in xrange(1, 23)])
                else:
                    ref_names = (["ref", "x", "X", "y", "Y", "m", "M"] +
                        [str(x) for x in xrange(1, 23)])
                RealTimeLogger.get().info("REF_NAMES: {}".format(ref_names))
                ref_options = []
                for name in ref_names:
                    # Put each in a -r option to retain the path
                    ref_options.append("-r")
                    ref_options.append(name)

                # Retain only the specified paths (only one should really exist)
                command = ['mod', '-N'] + ref_options + ['-t', str(job.cores), os.path.basename(graph_filename)]
                docker_call(work_dir=work_dir, parameters=command,
                            tool='quay.io/ucsc_cgl/vg:latest',
                            inputs=[graph_filename],
                            outfile=to_index_file)
                
        time.sleep(1)

        # Now we have the combined to-index graph in one vg file. We'll load
        # it (which deduplicates nodes/edges) and then find kmers.
        RealTimeLogger.get().info("Finding kmers in {} to {}".format(
            to_index_filename, kmers_filename))

        # Make the GCSA2 kmers file
        with open(kmers_filename, "w") as kmers_file:
            command = ['vg view -v {}'.format(os.path.basename(to_index_filename)),
                       'vg kmers -g -B -k {} -H 1000000000 -T 1000000001 -t {} -'.format(str(options.kmer_size), str(job.cores))]
            docker_call(work_dir=work_dir, parameters=command,
                        tools='quay.io/ucsc_cgl/vg:latest',
                        outfile=kmers_file)

        time.sleep(1)

        # Where do we put the GCSA2 index?
        gcsa_filename = graph_filename + ".gcsa"

        RealTimeLogger.get().info("GCSA-indexing {} to {}".format(
                kmers_filename, gcsa_filename))

        # Make the gcsa2 index. Make sure to use 3 doubling steps to work
        # around <https://github.com/vgteam/vg/issues/301>
        command = ['index', '-t', str(job.cores), '-i', 
            os.path.basename(kmers_filename), '-g', os.path.basename(gcsa_filename),
            '-X', '3', '-Z', '2000']
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest')

        # Where do we put the XG index?
        xg_filename = graph_filename + ".xg"

        RealTimeLogger.get().info("XG-indexing {} to {}".format(
                graph_filename, xg_filename))
        
        command = ['index', '-t', str(job.cores), '-x', 
            os.path.basename(xg_filename), os.path.basename(graph_filename)]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest')

    else:
        raise RuntimeError("Invalid indexing mode: " + options.index_mode)

    # Define a file to keep the compressed index in, so we can send it to
    # the output store.
    index_dir_tgz = "{}/index.tar.gz".format(
        job.fileStore.getLocalTempDir())

    # Now save the indexed graph directory to the file store. It can be
    # cleaned up since only our children use it.
    RealTimeLogger.get().info("Compressing index of {}".format(
        graph_filename))
    index_dir_id = write_global_directory(job.fileStore, graph_dir,
        cleanup=True, tee=index_dir_tgz)

    # Save it as output
    RealTimeLogger.get().info("Uploading index of {}".format(
        graph_filename))
    index_key = os.path.basename(index_dir_tgz)
    out_store.write_output_file(index_dir_tgz, index_key)
    RealTimeLogger.get().info("Index {} uploaded successfully".format(
        index_key))

    return index_key, index_dir_id

def run_only_indexing(job, options, inputGraphFileID):
    """ run chunked_call logic by itself.  todo-modularize existing pipeline enough
    so that we don't need to duplicate this code """

    RealTimeLogger.get().info("Uploading files to IO store")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)   
    out_store = IOStore.get(options.out_store)

    input_graph_basename = os.path.basename(options.vg_graph)
    RealTimeLogger.get().info("Uploading {} to {} on IO store".format(inputGraphFileID, input_graph_basename))
    fi = job.fileStore.readGlobalFile(inputGraphFileID)
    input_store.write_output_file(fi, input_graph_basename)

    index_key, index_dir_id = build_indexes(job, options)

    return index_key

def fetch_output_index(options, index_key):
    """ copy from out store to out dir 
    to do - having both outstore and out_dir seems redundant """
    
    # Create output directory if it doesn't exist
    try:
        os.makedirs(options.out_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: raise

    RealTimeLogger.get().info("Downloading {} to {}".format(index_key, options.out_dir))

    # Derive run_calling output files from its return value
    out_store = IOStore.get(options.out_store)

    # Read them out of the output store
    out_store.read_input_file(index_key, os.path.join(options.out_dir, index_key))


def main():
    """
    Wrapper for vg indexing. 
    """

    RealTimeLogger.start_master()
    options = parse_args() # This holds the nicely-parsed options object

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:
            
            # Upload local files to the remote IO Store
            inputGraphFileID = toil.importFile('file://'+options.vg_graph)
            
            # Make a root job
            root_job = Job.wrapJobFn(run_only_indexing, options, inputGraphFileID,
                                     cores=2, memory="5G", disk="2G")
            
            # Run the job and store the returned list of output files to download
            index_key = toil.start(root_job)
        else:
            index_key = toil.restart()
            
        # copy the indexes out of the output store and into options.out_dir
        fetch_output_index(options, index_key)

    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    RealTimeLogger.stop_master()

if __name__ == "__main__" :
    try:
        main()
    except Exception as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)

