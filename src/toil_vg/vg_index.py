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
from toil_vg.vg_common import *

def index_subparser(parser):
    """
    Create a subparser for indexing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("vg_graph", type=str,
        help="Input vg graph file path")
    parser.add_argument("out_store",
        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add indexing options
    index_parse_args(parser)

    # Add common docker options
    add_docker_tool_parse_args(parser)


def index_parse_args(parser):
    """ centralize indexing parameters here """
    
    parser.add_argument("--edge_max", type=int,
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int,
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--include_pruned", action="store_true",
        help="use the pruned graph in the index")
    parser.add_argument("--index_cores", type=int,
        help="number of threads during the indexing step")

def run_gcsa_indexing(job, options, inputGraphFileID):
    """
    Make the gcsa2 index. Return its store id
    """
    
    RealTimeLogger.get().info("Starting gcsa indexing...")
    start_time = timeit.default_timer()
    
    graph_name = os.path.basename(options.vg_graph)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Our local copy of the graph
    graph_filename = os.path.join(work_dir, graph_name)
    read_from_store(job, options, inputGraphFileID, graph_filename)    
    
    # Now run the indexer.
    RealTimeLogger.get().info("GCSA Indexing {}".format(options.vg_graph))
            
    # We want a GCSA2/xg index. We have to prune the graph ourselves.
    # See <https://github.com/vgteam/vg/issues/286>.

    # What will we use as our temp combined graph file (containing only
    # the bits of the graph we want to index, used for deduplication)?
    to_index_filename = os.path.join(work_dir, "to_index.vg")

    # Where will we save the kmers?
    kmers_filename = os.path.join(work_dir, "kmers.tsv")

    with open(to_index_filename, "w") as to_index_file:
        if options.include_pruned:
            RealTimeLogger.get().info("Pruning {} to {}".format(
                graph_filename, to_index_filename))

            # Prune out hard bits of the graph
            # and complex regions
            # and short disconnected chunks
            command = [['vg', 'mod', '-p', '-l', str(options.kmer_size), '-t', str(job.cores), '-e',
                        str(options.edge_max), os.path.basename(graph_filename)]]
            command.append(['vg', 'mod', '-S', '-l', str(options.kmer_size * 2), '-t', str(job.cores)])
            options.drunner.call(command, work_dir=work_dir, outfile=to_index_file)

        # Then append in the primary path. Since we don't knoiw what
        # "it's called, we retain "ref" and all the 19", "6", etc paths
        # "from 1KG.

        RealTimeLogger.get().info(
            "Adding primary path to {}".format(to_index_filename))
        # See
        # https://github.com/vgteam/vg/issues/318#issuecomment-215102199

        ref_options = []
        for name in options.path_name:
            # Put each in a -r option to retain the path
            ref_options.append("-r")
            ref_options.append(name)

        # Retain only the specified paths (only one should really exist)
        if len(options.path_name) > 0:
            command = ['vg', 'mod', '-N'] + ref_options + ['-t', str(job.cores), os.path.basename(graph_filename)]
            options.drunner.call(command, work_dir=work_dir, 
                                 inputs=[graph_filename],
                                 outfile=to_index_file)
                
    time.sleep(1)

    # Now we have the combined to-index graph in one vg file. We'll load
    # it (which deduplicates nodes/edges) and then find kmers.
    RealTimeLogger.get().info("Finding kmers in {} to {}".format(
        to_index_filename, kmers_filename))

    # Make the GCSA2 kmers file
    with open(kmers_filename, "w") as kmers_file:
        command = ['vg', 'kmers',  os.path.basename(to_index_filename), '-g', '-B', '-k',
                   str(options.kmer_size), '-H', '1000000000', '-T', '1000000001', '-t', str(job.cores)]
        options.drunner.call(command, work_dir=work_dir,
                             outfile=kmers_file)

    time.sleep(1)

    # Where do we put the GCSA2 index?
    gcsa_filename = graph_filename + ".gcsa"

    RealTimeLogger.get().info("GCSA-indexing {} to {}".format(
            kmers_filename, gcsa_filename))

    # Make the gcsa2 index. Make sure to use 3 doubling steps to work
    # around <https://github.com/vgteam/vg/issues/301>
    command = ['vg', 'index', '-t', str(job.cores), '-i', 
        os.path.basename(kmers_filename), '-g', os.path.basename(gcsa_filename),
        '-X', '3', '-Z', '2000']
    options.drunner.call(command, work_dir=work_dir)

    # Checkpoint index to output store
    if not options.force_outstore or options.tool == 'index':
        write_to_store(job, options, os.path.join(work_dir, gcsa_filename), use_out_store = True)
        write_to_store(job, options, os.path.join(work_dir, gcsa_filename) + ".lcp", use_out_store = True)

    # Not in standalone Mode, then we write it to the file store
    if options.tool != 'index':
        gcsa_file_id = write_to_store(job, options, os.path.join(work_dir, gcsa_filename))
        lcp_file_id = write_to_store(job, options, os.path.join(work_dir, gcsa_filename) + ".lcp")
        return gcsa_file_id, lcp_file_id

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealTimeLogger.get().info("Finished GCSA index. Process took {} seconds.".format(run_time))



def run_xg_indexing(job, options, inputGraphFileID):
    """ Make the xg index and return its store id
    """
    
    RealTimeLogger.get().info("Starting xg indexing...")
    start_time = timeit.default_timer()
    
    graph_name = os.path.basename(options.vg_graph)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Our local copy of the graph
    graph_filename = os.path.join(work_dir, graph_name)
    read_from_store(job, options, inputGraphFileID, graph_filename)    

    # Where do we put the XG index?
    xg_filename = graph_filename + ".xg"

    # Now run the indexer.
    RealTimeLogger.get().info("XG Indexing {}".format(options.vg_graph))            

    command = ['vg', 'index', '-t', str(job.cores), '-x', 
        os.path.basename(xg_filename), os.path.basename(graph_filename)]
    options.drunner.call(command, work_dir=work_dir)

    # Checkpoint index to output store
    if not options.force_outstore or options.tool == 'index':
        write_to_store(job, options, os.path.join(work_dir, xg_filename), use_out_store = True)

    # Not in standalone Mode, then we write it to the file store
    if options.tool != 'index':
        xg_file_id = write_to_store(job, options, os.path.join(work_dir, xg_filename))
        return xg_file_id

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealTimeLogger.get().info("Finished XG index. Process took {} seconds.".format(run_time))

    

def run_indexing(job, options, inputGraphFileID):
    """ run indexing logic by itself.  Return pair of idx for xg and gcsa output index files  
    """

    gcsa_and_lcp_ids = job.addChildJobFn(run_gcsa_indexing, options, inputGraphFileID,
                                      cores=options.index_cores, memory=options.index_mem, disk=options.index_disk).rv()
    xg_index_id = job.addChildJobFn(run_xg_indexing, options, inputGraphFileID,
                                      cores=options.index_cores, memory=options.index_mem, disk=options.index_disk).rv()

    return xg_index_id, gcsa_and_lcp_ids


def index_main(options):
    """
    Wrapper for vg indexing. 
    """
    
    RealTimeLogger.start_master()

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))

    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'index'

    # Throw error if something wrong with IOStore string
    IOStore.get(options.out_store)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:
            
            # Upload local files to the remote IO Store
            inputGraphFileID = import_to_store(toil, options, options.vg_graph)
            
            # Make a root job
            root_job = Job.wrapJobFn(run_indexing, options, inputGraphFileID,
                                     cores=options.misc_cores,
                                     memory=options.misc_mem,
                                     disk=options.misc_disk)
            
            # Run the job and store the returned list of output files to download
            index_key_and_id = toil.start(root_job)
        else:
            index_key_and_id = toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    RealTimeLogger.stop_master()

