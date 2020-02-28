#!/usr/bin/env python2.7
"""
vg_augment.py: augment a vg graph to include variation from a GAM alignment

"""

import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, struct, socket, threading
import string
import getpass
import pdb
import gzip
import logging

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_chunk import *

logger = logging.getLogger(__name__)

def augment_subparser(parser):
    """
    Create a subparser for augmenting.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--gam", type=make_url, required=True,
                        help="GAM to augment")
    parser.add_argument("--graph", type=make_url, required=True,
                        help="graph to augment")
        
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add augment options
    augment_parse_args(parser, stand_alone = True)

    # Add common chunking options shared with vg_chunk
    chunk_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)


def augment_parse_args(parser, stand_alone = False):
    """
    Define map arguments shared with mapeval and run
    """
    if stand_alone:
        parser.add_argument("--augment_gam", action="store_true",
                            help="produce and augmented GAM")    
    parser.add_argument("--min_augment_coverage", type=int, 
                        help="minimum coverage for breakpoint to be applied")
    parser.add_argument("--expected_coverage", type=int,
                        help="expected coverage.  only affects memory usage.  use if coverage >> 100")
    parser.add_argument("--min_mapq", type=int,
                        help="ignore reads with MAPQ less than this")
    parser.add_argument("--min_baseq", type=int,
                        help="ignore edits with minimum average base quality less than this")
    parser.add_argument("--augment_cores", type=int,
                        help="number of threads during augmentation")
    parser.add_argument("--augment_mem", type=str,
                        help="memory alotment during augmentation")

def run_chunked_augmenting(job, context,
                           graph_id,
                           graph_basename,
                           gam_id,
                           gam_basename,
                           batch_input=None,
                           all_path_components=False,
                           chunk_paths=[],
                           connected_component_chunking=False,
                           output_format=None,
                           augment_gam=False,
                           min_augment_coverage=None,
                           expected_coverage=None,
                           min_mapq=None, 
                           min_baseq=None,
                           to_outstore=False):
    """
    Run a chunking job (if desired), then augment the results
    """

    # base case: only one input
    if batch_input is None:
        # chunk if necessary
        if all_path_components or connected_component_chunking or len(chunk_paths) > 1:
            child_job = Job()
            job.addChild(child_job)
            chunk_job = child_job.addChildJobFn(run_chunking, context,
                                                graph_id=graph_id,
                                                graph_basename=graph_basename,
                                                chunk_paths=chunk_paths,
                                                connected_component_chunking=connected_component_chunking,
                                                output_format=output_format,
                                                gam_id=gam_id,
                                                to_outstore=False,
                                                cores=context.config.chunk_cores,
                                                memory=context.config.chunk_mem,
                                                disk=context.config.chunk_disk)
            batch_input = chunk_job.rv()

            # recurse on chunks
            recurse_job = child_job.addFollowOnJobFn(run_chunked_augmenting, context,
                                                     graph_id=None,
                                                     graph_basename=None,
                                                     gam_id=None,
                                                     gam_basename=None,
                                                     batch_input=batch_input,
                                                     all_path_components=all_path_components,
                                                     chunk_paths=chunk_paths,
                                                     connected_component_chunking=connected_component_chunking,
                                                     output_format=output_format,
                                                     augment_gam=augment_gam,
                                                     min_augment_coverage=min_augment_coverage,
                                                     expected_coverage=expected_coverage,
                                                     min_mapq=min_mapq, 
                                                     min_baseq=min_baseq,
                                                     to_outstore=to_outstore)
            return recurse_job.rv()
        else:
            #phony up chunk output for single input
            batch_input = { 'all' : [graph_id, graph_basename] }
            if gam_id:
                batch_input['all'] += [gam_id, gam_basename]

    # run the augmenting on each chunk
    assert batch_input

    augment_results = []
    for chunk_name, chunk_results in list(batch_input.items()):
        augment_job = job.addChildJobFn(run_augmenting, context,
                                        graph_id=chunk_results[0],
                                        graph_basename=chunk_results[1],
                                        gam_id=chunk_results[2],
                                        gam_basename=chunk_results[3],
                                        augment_gam=augment_gam,
                                        min_augment_coverage=min_augment_coverage,
                                        expected_coverage=expected_coverage,
                                        min_mapq=min_mapq,
                                        min_baseq=min_baseq,
                                        to_outstore=to_outstore,
                                        cores=context.config.augment_cores,
                                        memory=context.config.augment_mem,
                                        disk=context.config.augment_disk)
        
        augment_results.append((chunk_name, augment_job.rv()))

    return augment_results
        

def run_augmenting(job, context,
                  graph_id,
                  graph_basename,
                  gam_id,
                  gam_basename,
                  augment_gam=False,
                  min_augment_coverage=None,
                  expected_coverage=None,
                  min_mapq=None,
                  min_baseq=None,
                  to_outstore=False):
    """
    Augment the graph (and gam if wanted)
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    graph_path = os.path.join(work_dir, graph_basename)
    job.fileStore.readGlobalFile(graph_id, graph_path)    
    gam_path = os.path.join(work_dir, gam_basename)
    job.fileStore.readGlobalFile(gam_id, gam_path)

    augment_cmd = ['vg', 'augment', os.path.basename(graph_path), os.path.basename(gam_path), '-t', str(job.cores)]

    # hardcoded naming convention: tack on an -aug to input paths
    graph_name, graph_ext = os.path.splitext(graph_basename)
    augmented_graph_path = os.path.join(work_dir, graph_name + '-aug' + graph_ext)
    augmented_gam_path = os.path.join(work_dir, remove_ext(gam_basename, '.gam') + '-aug.gam')
    if augment_gam:
        augment_cmd += ['-A', os.path.basename(augmented_gam_path)]

    # optional stuff
    if min_augment_coverage is not None:
        augment_cmd += ['-m', str(min_augment_coverage)]
    if expected_coverage is not None:
        augment_cmd += ['-c', str(expected_coverage)]
    if min_mapq is not None:
        augment_cmd += ['-Q', str(min_mapq)]
    if min_baseq is not None:
        augment_cmd += ['-q', str(min_baseq)]
    if context.config.augment_opts:
        augment_cmd += context.config.augment_opts
    # always support subgraphs
    augment_cmd += ['-s']

    # run the command
    try:
        with open(augmented_graph_path, 'wb') as augmented_graph_file:
            context.runner.call(job, augment_cmd, work_dir = work_dir, outfile = augmented_graph_file)
    except Exception as e:
        logging.error("Augment failed. Dumping input files to outstore.")
        for dump_path in [graph_path, gam_path]:
            if dump_path and os.path.isfile(dump_path):
                context.write_output_file(job, dump_path)
        raise

    # return the output
    write_fn = context.write_output_file if to_outstore else context.write_intermediate_file
    aug_graph_id = write_fn(job, augmented_graph_path)
    aug_gam_id = write_fn(job, augmented_gam_path) if augment_gam else None
    return aug_graph_id, aug_gam_id
    
def augment_main(context, options):
    """
    Wrapper for vg augment.
    """

    validate_chunk_options(options, chunk_optional=True)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)

            # Upload local files to the job store
            inputGraphFileID = importer.load(options.graph)
            inputGamFileID = importer.load(options.gam)

            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_chunked_augmenting, context,
                                     graph_id = importer.resolve(inputGraphFileID),
                                     graph_basename = os.path.basename(options.graph),
                                     gam_id = importer.resolve(inputGamFileID),
                                     gam_basename = os.path.basename(options.gam),
                                     batch_input=None,
                                     all_path_components=options.ref_path_chunking,
                                     chunk_paths=options.ref_paths,
                                     connected_component_chunking=options.connected_component_chunking,
                                     output_format=options.output_format,
                                     augment_gam=options.augment_gam,
                                     min_augment_coverage=options.min_augment_coverage,
                                     expected_coverage=options.expected_coverage,
                                     min_mapq=options.min_mapq,
                                     min_baseq=options.min_baseq,
                                     to_outstore=True)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
            init_job.addFollowOn(root_job)            
            
            # Run the job and store the returned list of output files to download
            toil.start(init_job)
        else:
            toil.restart()

    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
