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
from toil_lib import require

def index_subparser(parser):
    """
    Create a subparser for indexing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
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

    parser.add_argument("--graphs", nargs='+', required=True,
                        help="input graph(s). one per chromosome (separated by space)")

    parser.add_argument("--index_cores", type=int,
        help="number of threads during the indexing step")

def compute_gcsa_kmers(options, cores, work_dir, graph_filename, chrom):
    """
    Make the kmers file, return the path
    """
    RealTimeLogger.get().info("Starting gcsa kmers...")
    start_time = timeit.default_timer()

    # Output
    output_kmers_filename = graph_filename + '.kmers'
    
    # Place where we put pruned vg
    to_index_filename = graph_filename
    
    if len(options.prune_opts) > 0:
        to_index_filename = os.path.join(work_dir, "to_index.vg")
        with open(to_index_filename, "w") as to_index_file:
            pcores = cores if len(options.prune_opts_2) == 0 else cores / 2
            command = ['vg', 'mod', os.path.basename(graph_filename), '-t',
                       str(pcores)] + options.prune_opts
            # tack on 2nd vg mod command if specified
            # note: perhaps this is a bakeoff relic that we don't need anymore
            # if we push -S to first command. 
            if len(options.prune_opts_2) > 0:
                command = [command]
                command.append(['vg', 'mod', '-', '-t', str(pcores)] + options.prune_opts_2)
            options.drunner.call(command, work_dir=work_dir, outfile=to_index_file)

            # Then append in the primary path.
            command = ['vg', 'mod', '-N', '-r', chrom]
            command += ['-t', str(cores), os.path.basename(graph_filename)]
            options.drunner.call(command, work_dir=work_dir, outfile=to_index_file)

    # Now we have the combined to-index graph in one vg file. We'll load
    # it (which deduplicates nodes/edges) and then find kmers.
    RealTimeLogger.get().info("Finding kmers in {} to {}".format(
        to_index_filename, output_kmers_filename))

    # Make the GCSA2 kmers file
    with open(output_kmers_filename, "w") as kmers_file:
        command = ['vg', 'kmers',  os.path.basename(to_index_filename), '-t', str(cores)]
        command += options.kmers_opts
        options.drunner.call(command, work_dir=work_dir, outfile=kmers_file)

    # Dont' need to keep pruned graph around after kmers computed.
    if to_index_filename != graph_filename:
        os.remove(to_index_filename)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealTimeLogger.get().info("Finished GCSA kmers. Process took {} seconds.".format(run_time))

    return output_kmers_filename
    
def run_gcsa_indexing(job, options, input_graph_ids):
    """
    Make the gcsa2 index. Return its store id
    """
    
    RealTimeLogger.get().info("Starting gcsa indexing...")
    start_time = timeit.default_timer()     

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    kmers_filenames = []
    
    # Compute our kmers for each input graph (in series)
    # is it worth it to distrbute?  files are so big to move around...
    for graph_i, input_graph_id in enumerate(input_graph_ids):
        graph_filename = os.path.join(work_dir, os.path.basename(options.graphs[graph_i]))
        read_from_store(job, options, input_graph_id, graph_filename)
        kmers_filenames.append(compute_gcsa_kmers(options, job.cores, work_dir,
                                                  graph_filename, options.chroms[graph_i]))

    # Where do we put the GCSA2 index?
    gcsa_filename = "genome.gcsa"

    command = ['vg', 'index', '-g', os.path.basename(gcsa_filename)] + options.gcsa_opts
    command += ['-t', str(job.cores)]
    for kmers_filename in kmers_filenames:
        command += ['-i', os.path.basename(kmers_filename)]
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

def run_xg_indexing(job, options, inputGraphFileIDs):
    """ Make the xg index and return its store id
    """
    
    RealTimeLogger.get().info("Starting xg indexing...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Our local copy of the graphs
    graph_filenames = []
    for i, graph_id in enumerate(inputGraphFileIDs):
        graph_filename = os.path.join(work_dir, '{}.vg'.format(options.chroms[i]))
        read_from_store(job, options, graph_id, graph_filename)
        graph_filenames.append(os.path.basename(graph_filename))

    # Where do we put the XG index?
    xg_filename = graph_filename + ".xg"

    # Now run the indexer.
    RealTimeLogger.get().info("XG Indexing {}".format(str(graph_filenames)))            

    command = ['vg', 'index', '-t', str(job.cores), '-x', os.path.basename(xg_filename)]
    command += graph_filenames
    
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

    

def run_indexing(job, options, inputGraphFileIDs):
    """ run indexing logic by itself.  Return pair of idx for xg and gcsa output index files  
    """

    gcsa_and_lcp_ids = job.addChildJobFn(run_gcsa_indexing, options, inputGraphFileIDs,
                                      cores=options.index_cores, memory=options.index_mem, disk=options.index_disk).rv()
    xg_index_id = job.addChildJobFn(run_xg_indexing, options, inputGraphFileIDs,
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

    require(len(options.chroms) == len(options.graphs), '--chrom and --graph must have'
            ' same number of arguments')

    # Throw error if something wrong with IOStore string
    IOStore.get(options.out_store)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:
            
            # Upload local files to the remote IO Store
            inputGraphFileIDs = []
            for graph in options.graphs:
                inputGraphFileIDs.append(import_to_store(toil, options, graph))
            
            # Make a root job
            root_job = Job.wrapJobFn(run_indexing, options, inputGraphFileIDs,
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

