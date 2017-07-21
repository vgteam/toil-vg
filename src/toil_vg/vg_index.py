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
import logging

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context

logger = logging.getLogger(__name__)

def index_subparser(parser):
    """
    Create a subparser for indexing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")

    parser.add_argument("--skip_xg", action="store_true",
                        help="Do not generate xg index")
    parser.add_argument("--skip_gcsa", action="store_true",
                        help="Do not generate gcsa index")
    parser.add_argument("--skip_id_ranges", action="store_true",
                        help="Do not generate id_ranges.tsv")
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add indexing options
    index_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)


def index_parse_args(parser):
    """ centralize indexing parameters here """

    parser.add_argument("--graphs", nargs='+', type=make_url,
                        help="input graph(s). one per chromosome (separated by space)")

    parser.add_argument("--chroms", nargs='+',
                        help="Name(s) of reference path in graph(s) (separated by space).  If --graphs "
                        " specified, must be same length/order as --chroms")

    parser.add_argument("--gcsa_index_cores", type=int,
        help="number of threads during the gcsa indexing step")

    parser.add_argument("--kmers_cores", type=int,
        help="number of threads during the gcsa kmers step")

    parser.add_argument("--index_name", type=str, default='index',
                        help="name of index files. <name>.xg, <name>.gcsa etc.")

    parser.add_argument("--prune_opts", type=str,
                        help="Options to pass to vg mod for pruning phase.")
    parser.add_argument("--kmers_opts", type=str,
                        help="Options to pass to vg kmers.")
    parser.add_argument("--gcsa_opts", type=str,
                        help="Options to pass to gcsa indexing.")

    parser.add_argument("--vcf_phasing", type=make_url,
                        help="Import phasing information from VCF into xg")

def validate_index_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """                           
    require(options.chroms and options.graphs, '--chroms and --graphs must be specified')
    require(len(options.chroms) == len(options.graphs), '--chroms and --graphs must have'
            ' same number of arguments')
    if options.vcf_phasing:
        require(options.vcf_phasing.endswith('.vcf.gz'), 'input phasing file must end with .vcf.gz')
    
def run_gcsa_prune(job, context, graph_name, input_graph_id, primary_paths=[]):
    """
    Make the pruned graph, do kmers as a follow up and return kmers id.
    Retains the specified primary paths.
    """
    RealtimeLogger.info("Starting graph-pruning for gcsa kmers...")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Place where we put pruned vg
    to_index_filename = None
    
    if len(context.config.prune_opts) > 0:

        # Download input graph
        graph_filename = os.path.join(work_dir, graph_name)
        job.fileStore.readGlobalFile(input_graph_id, graph_filename)

        to_index_filename = os.path.join(work_dir, "pruned_{}".format(graph_name))
        with open(to_index_filename, "w") as to_index_file:
            command = ['vg', 'mod', os.path.basename(graph_filename), '-t',
                       str(job.cores)] + context.config.prune_opts
            # tack on 2nd vg mod command if specified
            # note: perhaps this is a bakeoff relic that we don't need anymore
            # if we push -S to first command. 
            if len(context.config.prune_opts_2) > 0:
                command = [command]
                command.append(['vg', 'mod', '-', '-t', str(job.cores)] + context.config.prune_opts_2)
            context.runner.call(job, command, work_dir=work_dir, outfile=to_index_file)
            
            # Then append in the primary path.
            command = ['vg', 'mod', '-N']
            for primary_path in primary_paths:
                # Send along all the primary paths
                command.append('-r')
                command.append(primary_path)
            command += ['-t', str(job.cores), os.path.basename(graph_filename)]
            context.runner.call(job, command, work_dir=work_dir, outfile=to_index_file)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished pruning. Process took {} seconds.".format(run_time))

    if to_index_filename is None:
        # no pruning done: just pass along the input graph as is
        pruned_graph_id = input_graph_id
    else:
        pruned_graph_id = context.write_intermediate_file(job, to_index_filename)
    
    return job.addFollowOnJobFn(run_gcsa_kmers, context, graph_name,
                                pruned_graph_id, 
                                cores=context.config.kmers_cores, memory=context.config.kmers_mem,
                                disk=context.config.kmers_disk).rv()

def run_gcsa_kmers(job, context, graph_name, input_graph_id):
    """
    Make the kmers file, return its id
    """
    RealtimeLogger.info("Starting gcsa kmers...")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download input graph
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(input_graph_id, graph_filename)

    # Output
    output_kmers_filename = graph_filename + '.kmers'
   
    RealtimeLogger.info("Finding kmers in {} to {}".format(
        graph_filename, output_kmers_filename))

    # Make the GCSA2 kmers file
    with open(output_kmers_filename, "w") as kmers_file:
        command = ['vg', 'kmers',  os.path.basename(graph_filename), '-t', str(job.cores)]
        command += context.config.kmers_opts
        context.runner.call(job, command, work_dir=work_dir, outfile=kmers_file)

    # Back to store
    output_kmers_id = context.write_intermediate_file(job, output_kmers_filename)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GCSA kmers. Process took {} seconds.".format(run_time))

    return output_kmers_id

def run_gcsa_prep(job, context, input_graph_ids,
                  graph_names, index_name, chroms,
                  primary_path_override=None):
    """
    Do all the preprocessing for gcsa indexing (pruning and kmers)
    Then launch the indexing as follow-on
    """    
    RealtimeLogger.info("Starting gcsa preprocessing...")
    start_time = timeit.default_timer()     

    kmers_ids = []
    
    # Compute our kmers for each input graph (in series)
    # is it worth it to distrbute?  files are so big to move around...
    for graph_i, input_graph_id in enumerate(input_graph_ids):
        # For each input graph
        
        # Determine the primary path list to use
        primary_paths = ([chroms[graph_i]] if primary_path_override
            is None else primary_path_override)
        
        # Make the kmers, passing along the primary path names
        kmers_id = job.addChildJobFn(run_gcsa_prune, context, graph_names[graph_i], input_graph_id,
                                     primary_paths=primary_paths,
                                     cores=context.config.prune_cores, memory=context.config.prune_mem,
                                     disk=context.config.prune_disk).rv()
        kmers_ids.append(kmers_id)

    return job.addFollowOnJobFn(run_gcsa_indexing, context, kmers_ids,
                                graph_names, index_name,
                                cores=context.config.gcsa_index_cores,
                                memory=context.config.gcsa_index_mem,
                                disk=context.config.gcsa_index_disk).rv()
    
def run_gcsa_indexing(job, context, kmers_ids, graph_names, index_name):
    """
    Make the gcsa2 index. Return its store id
    """
    
    RealtimeLogger.info("Starting gcsa indexing...")
    start_time = timeit.default_timer()     

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download all the kmers.  
    kmers_filenames = []
    
    for graph_i, kmers_id in enumerate(kmers_ids):
        kmers_filename = os.path.join(work_dir, os.path.basename(graph_names[graph_i]) + '.kmers')
        job.fileStore.readGlobalFile(kmers_id, kmers_filename)
        kmers_filenames.append(kmers_filename)

    # Where do we put the GCSA2 index?
    gcsa_filename = "{}.gcsa".format(index_name)

    command = ['vg', 'index', '-g', os.path.basename(gcsa_filename)] + context.config.gcsa_opts
    command += ['-t', str(job.cores)]
    for kmers_filename in kmers_filenames:
        command += ['-i', os.path.basename(kmers_filename)]
    context.runner.call(job, command, work_dir=work_dir)

    # Checkpoint index to output store
    gcsa_file_id = context.write_output_file(job, os.path.join(work_dir, gcsa_filename))
    lcp_file_id = context.write_output_file(job, os.path.join(work_dir, gcsa_filename) + ".lcp")

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GCSA index. Process took {} seconds.".format(run_time))

    return gcsa_file_id, lcp_file_id


def run_xg_indexing(job, context, inputGraphFileIDs, graph_names, index_name,
                    vcf_phasing_file_id = None, tbi_phasing_file_id = None):
    """ Make the xg index and return its store id
    """
    
    RealtimeLogger.info("Starting xg indexing...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Our local copy of the graphs
    graph_filenames = []
    for i, graph_id in enumerate(inputGraphFileIDs):
        graph_filename = os.path.join(work_dir, graph_names[i])
        job.fileStore.readGlobalFile(graph_id, graph_filename)
        graph_filenames.append(os.path.basename(graph_filename))

    # Get the vcf file for making gpbwt
    if vcf_phasing_file_id:
        phasing_file = os.path.join(work_dir, 'phasing.vcf.gz')
        job.fileStore.readGlobalFile(vcf_phasing_file_id, phasing_file)
        job.fileStore.readGlobalFile(tbi_phasing_file_id, phasing_file + '.tbi')
        phasing_opts = ['-v', os.path.basename(phasing_file)]
    else:
        phasing_opts = []

    # Where do we put the XG index?
    xg_filename = "{}.xg".format(index_name)

    # Now run the indexer.
    RealtimeLogger.info("XG Indexing {}".format(str(graph_filenames)))

    command = ['vg', 'index', '-t', str(job.cores), '-x', os.path.basename(xg_filename)]
    command += phasing_opts + graph_filenames
    
    context.runner.call(job, command, work_dir=work_dir)

    # Checkpoint index to output store
    xg_file_id = context.write_output_file(job, os.path.join(work_dir, xg_filename))

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished XG index. Process took {} seconds.".format(run_time))

    return xg_file_id


def run_id_ranges(job, context, inputGraphFileIDs, graph_names, index_name, chroms):
    """ Make a file of chrom_name <tab> first_id <tab> last_id covering the 
    id ranges of all chromosomes.  This is to speed up gam splitting down the road. 
    """
    
    RealtimeLogger.info("Starting id ranges...")
    start_time = timeit.default_timer()
    
    # Our id ranges (list of triples)
    id_ranges = []

    # Get the range for one graph per job. 
    for graph_id, graph_name, chrom in zip(inputGraphFileIDs, graph_names, chroms):
        id_range = job.addChildJobFn(run_id_range, context, graph_id, graph_name, chrom,
                                     cores=context.config.prune_cores,
                                     memory=context.config.prune_mem, disk=context.config.prune_disk).rv()
        
        id_ranges.append(id_range)

    # Merge them into a file and return its id
    return job.addFollowOnJobFn(run_merge_id_ranges, context, id_ranges, index_name,
                                cores=context.config.misc_cores, memory=context.config.misc_mem,
                                disk=context.config.misc_disk).rv()

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished id ranges. Process took {} seconds.".format(run_time))
    
def run_id_range(job, context, graph_id, graph_name, chrom):
    """
    Compute a node id range for a graph (which should be an entire contig/chromosome with
    contiguous id space -- see vg ids) using vg stats
    """
    work_dir = job.fileStore.getLocalTempDir()

    # download graph
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(graph_id, graph_filename)

    #run vg stats
    #expect result of form node-id-range <tab> first:last
    command = ['vg', 'stats', '-r', os.path.basename(graph_filename)]
    stats_out = context.runner.call(job, command, work_dir=work_dir, check_output = True).strip().split()
    assert stats_out[0] == 'node-id-range'
    first, last = stats_out[1].split(':')

    return chrom, first, last
    
def run_merge_id_ranges(job, context, id_ranges, index_name):
    """ create a BED-style file of id ranges
    """
    work_dir = job.fileStore.getLocalTempDir()

    # Where do we put the XG index?
    id_range_filename = os.path.join(work_dir, '{}_id_ranges.tsv'.format(index_name))

    with open(id_range_filename, 'w') as f:
        for id_range in id_ranges:
            f.write('{}\t{}\t{}\n'.format(*id_range))

    # Checkpoint index to output store
    return context.write_output_file(job, id_range_filename)

def run_indexing(job, context, inputGraphFileIDs,
                 graph_names, index_name, chroms,
                 vcf_phasing_file_id = None, tbi_phasing_file_id = None,
                 skip_xg=False, skip_gcsa=False, skip_id_ranges=False):
    """ run indexing logic by itself.  Return pair of idx for xg and gcsa output index files  
    """

    if not skip_gcsa:
        gcsa_and_lcp_ids = job.addChildJobFn(run_gcsa_prep, context, inputGraphFileIDs,
                                             graph_names, index_name, chroms,
                                             cores=context.config.misc_cores,
                                             memory=context.config.misc_mem,
                                             disk=context.config.misc_disk).rv()
    else:
        gcsa_and_lcp_ids = None
    if not skip_xg:
        xg_index_id = job.addChildJobFn(run_xg_indexing, context, inputGraphFileIDs,
                                        graph_names, index_name,
                                        vcf_phasing_file_id, tbi_phasing_file_id,
                                        cores=context.config.xg_index_cores,
                                        memory=context.config.xg_index_mem,
                                        disk=context.config.xg_index_disk).rv()
    else:
        xg_index_id = None
        
    if len(inputGraphFileIDs) > 1 and not skip_id_ranges:
        id_ranges_id = job.addChildJobFn(run_id_ranges, context, inputGraphFileIDs,
                                         graph_names, index_name, chroms,
                                         cores=context.config.misc_cores,
                                         memory=context.config.misc_mem,
                                         disk=context.config.misc_disk).rv()
    else:
        id_ranges_id = None

    return xg_index_id, gcsa_and_lcp_ids, id_ranges_id


def index_main(context, options):
    """
    Wrapper for vg indexing. 
    """

    # check some options
    validate_index_options(options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            inputGraphFileIDs = []
            for graph in options.graphs:
                inputGraphFileIDs.append(toil.importFile(graph))
            if options.vcf_phasing:
                inputPhasingVCFFileID = toil.importFile(options.vcf_phasing)
                inputPhasingTBIFileID = toil.importFile(options.vcf_phasing + '.tbi')
            else:
                inputPhasingVCFFileID = None
                inputPhasingTBIFileID = None

            # Handy to have meaningful filenames throughout, so we remember
            # the input graph names
            graph_names = [os.path.basename(i) for i in options.graphs]

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_indexing, context, inputGraphFileIDs,
                                     graph_names, options.index_name, options.chroms,
                                     inputPhasingVCFFileID, inputPhasingTBIFileID,
                                     options.skip_xg, options.skip_gcsa, options.skip_id_ranges,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
            
            # Run the job and store the returned list of output files to download
            index_key_and_id = toil.start(root_job)
        else:
            index_key_and_id = toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
