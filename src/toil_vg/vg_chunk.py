#!/usr/bin/env python2.7
"""
vg_chunk.py: split a graph and/or GAM into chunks by connected component

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

logger = logging.getLogger(__name__)

def chunk_subparser(parser):
    """
    Create a subparser for augmenting.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options    
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--gam", type=make_url,
                        help="GAM to chunk")
    parser.add_argument("--graph", type=make_url, required=True,
                        help="graph or xg index to chunk")
        
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add augment options
    chunk_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)


def chunk_parse_args(parser, path_components=True):
    """
    Define chunk arguments that may be shared with other commands
    """

    parser.add_argument("--connected_component_chunking", action="store_true",
                        help="split into connected components")
    parser.add_argument("--ref_path_chunking", action="store_true",
                        help="chunk on --ref_paths if specified, all paths otherwise")
    parser.add_argument("--ref_paths", nargs='+', default=[],
                        help="reference paths to call (and chunk) on")    
    parser.add_argument("--output_format", choices=["pg", "hg", "vg"], default="pg",
                        help="output format [pg]")
    parser.add_argument("--chunk_cores", type=int,
                        help="number of threads used for extracting chunks for calling")
    parser.add_argument("--chunk_mem", type=str,
                        help="memory alotment for extracting chunks for calling")

    
def validate_chunk_options(options, chunk_optional=False):
    num_opts = [options.connected_component_chunking, options.ref_path_chunking].count(True)
    if chunk_optional == False:
        require(num_opts == 1,
                "Must specify (exactly) one of --connected_component_chunking, ref_path_chunking")
        if options.ref_paths:
            require(options.ref_path_chunking == True,
                "Must specify --ref_path_chunking when using --ref_paths")
    else:
        require(num_opts in [0, 1],
                "Must specify at most one of --connected_component_chunking or --ref_path_chunking")

    if options.gam:
        require(options.gam.endswith('.gam'),
                "Input GAM file must have .gam extension")

def run_chunking(job, context,
                 graph_id,
                 graph_basename,
                 chunk_paths,
                 connected_component_chunking,
                 output_format,
                 gam_id = None,
                 to_outstore = False):
    """ thin wrapper of vg chunk.  it will return a map from path name to id of chunked file"""

    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    graph_path = os.path.join(work_dir, graph_basename)
    job.fileStore.readGlobalFile(graph_id, graph_path)
    input_opts = ['-x', os.path.basename(graph_path)]
    if gam_id:
        gam_path = os.path.join(work_dir, 'aln.gam')
        job.fileStore.readGlobalFile(gam_id, gam_path)
        input_opts += ['-a', os.path.basename(gam_path), '-g']

    paths_path = os.path.join(work_dir, 'paths.txt')        
    if chunk_paths:
        with open(paths_path, 'w') as path_file:
            for chunk_path in chunk_paths:
                path_file.write(chunk_path + '\n')
        input_opts += ['-P', os.path.basename(paths_path)]

    # output options
    chunk_prefix = 'chunk/{}'.format(os.path.splitext(graph_basename)[0])
    os.makedirs(os.path.join(work_dir, os.path.dirname(chunk_prefix)))    
    output_opts = ['-b', chunk_prefix, ]
    output_bed_path = os.path.join(work_dir, 'chunks.bed')
    output_opts += ['-E', os.path.basename(output_bed_path)]
    output_opts += ['-O', output_format]

    # general options
    if connected_component_chunking or len(chunk_paths) > 0:
        gen_opts = ['-C']
    else:
        gen_opts = ['-M']
    gen_opts += ['-t', str(job.cores)]

    # Run vg chunk
    try:
        context.runner.call(job, ['vg', 'chunk'] + gen_opts + input_opts + output_opts,
                            work_dir = work_dir)
    except Exception as e:
        logging.error("Chunk failed. Dumping input files to outstore.")
        for dump_path in [graph_path, gam_path, paths_path]:
            if dump_path and os.path.isfile(dump_path):
                context.write_output_file(job, dump_path)
        raise
        
    # Scrape the BED into dictionary that maps path name to file id
    chunk_output = {}
    write_fn = context.write_output_file if to_outstore else context.write_intermediate_file
    with open(output_bed_path) as output_bed:
        for line in output_bed:
            toks = line.split('\t')
            if len(toks) > 3:
                graph_chunk_path = os.path.join(work_dir, toks[3].rstrip())
                # be robust to the vagaries of vg chunk: deal with graph or gam extension in bed
                if graph_chunk_path.endswith('.gam'):
                    graph_chunk_path = remove_ext(graph_chunk_path, '.gam') + '.' + output_format
                graph_chunk_id = write_fn(job, graph_chunk_path)
                chunk_output[toks[0]] = [graph_chunk_id, os.path.basename(graph_chunk_path)]
                if gam_id:
                    gam_chunk_path = os.path.join(work_dir, toks[3].rstrip())
                    if not gam_chunk_path.endswith('.gam'):
                        gam_chunk_path = remove_ext(gam_chunk_path, '.' + output_format) + '.gam'
                    # vg chunk's not going to write empty files, so make sure we have one
                    open(gam_chunk_path, 'a').close()
                    gam_chunk_id = write_fn(job, gam_chunk_path)
                    chunk_output[toks[0]] += [gam_chunk_id, os.path.basename(gam_chunk_path)]

    return chunk_output
    
    
def chunk_main(context, options):
    """ entrypoint for calling """

    validate_chunk_options(options)
            
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)

            # Upload local files to the job store
            inputGraphFileID = importer.load(options.graph)
            inputGamFileID = None
            if options.gam:
                inputGamFileID = importer.load(options.gam)

            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_chunking, context,
                                     importer.resolve(inputGraphFileID),
                                     os.path.basename(options.graph),
                                     chunk_paths=options.ref_paths,
                                     connected_component_chunking=options.connected_component_chunking,
                                     output_format=options.output_format,
                                     gam_id = importer.resolve(inputGamFileID),
                                     to_outstore = True,
                                     cores=context.config.chunk_cores,
                                     memory=context.config.chunk_mem,
                                     disk=context.config.chunk_disk)

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
    
    
