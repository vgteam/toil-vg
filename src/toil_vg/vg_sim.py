#!/usr/bin/env python2.7
"""
vg_sim.py: this wrapper to run vg sim in parallel

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
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

logger = logging.getLogger(__name__)

def sim_subparser(parser):
    """
    Create a subparser for mapping.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("xg_index", type=str,
                        help="Path to xg index")
    parser.add_argument("num_reads", type=int,
                        help="Number of reads to simulate")
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--gam", action="store_true",
                        help="Output GAM file, annotated gam file and truth positions")
    parser.add_argument("--sim_opts", type=str,
                        help="arguments for vg sim (wrapped in \"\"). Do not include -x, -n, -s, or -a")
    parser.add_argument("--sim_chunks", type=int, default=1,
                        help="split simulation into this many chunks, run in parallel when possible")
    parser.add_argument("--seed", type=int, default=None,
                        help="random seed")
                        
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def run_sim(job, options, xg_file_id):
    """  
    run a bunch of simulation child jobs, merge up their output as a follow on
    """

    # each element is either reads_chunk_id or (gam_chunk_id, annot_gam_chunk_id, true_pos_chunk_id)
    # if --gam not specified
    sim_out_id_infos = []
    for chunk_i in range(options.sim_chunks):
        num_reads = options.num_reads / options.sim_chunks
        if chunk_i == options.sim_chunks - 1:
            num_reads += options.num_reads % options.sim_chunks
        sim_out_id_info = job.addChildJobFn(run_sim_chunk, options,xg_file_id,  chunk_i, num_reads,
                                            cores=options.sim_cores, memory=options.sim_mem,
                                            disk=options.sim_disk).rv()
        sim_out_id_infos.append(sim_out_id_info)

    return job.addFollowOnJobFn(run_merge_sim_chunks, options, sim_out_id_infos,
                                cores=options.sim_cores, memory=options.sim_mem,
                                disk=options.sim_disk).rv()

def run_sim_chunk(job, options, xg_file_id, chunk_i, num_reads):
    """
    simulate some reads (and optionally gam),
    return either reads_chunk_id or (gam_chunk_id, annot_gam_chunk_id, true_pos_chunk_id)
    if --gam specified
    """

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # read the xg file
    xg_file = os.path.join(work_dir, os.path.basename(options.xg_index))
    read_from_store(job, options, xg_file_id, xg_file)

    # run vg sim
    sim_cmd = ['vg', 'sim', '-x', os.path.basename(xg_file), '-n', num_reads] + options.sim_opts
    if options.seed is not None:
        sim_cmd += ['-s', options.seed + chunk_i]

    if not options.gam:
        # output reads
        reads_file = os.path.join(work_dir, 'sim_reads_{}'.format(chunk_i))

        # run vg sim
        with open(reads_file, 'w') as output_reads:
            options.drunner.call(job, sim_cmd, work_dir = work_dir, outfile=output_reads)

        # write to the store
        reads_chunk_id = write_to_store(job, options, reads_file)
        return reads_chunk_id
    else:
        # output gam
        gam_file = os.path.join(work_dir, 'sim_{}.gam'.format(chunk_i))
        gam_annot_file = os.path.join(work_dir, 'sim_{}_annot.gam'.format(chunk_i))
        gam_annot_json = os.path.join(work_dir, 'sim_{}_annot.json'.format(chunk_i))

        # run vg sim, write output gam, annotated gam, annotaged gam json
        # (from vg/scripts/map-sim)
        cmd = [sim_cmd + ['-a']]
        cmd.append(['tee', os.path.basename(gam_file)])
        cmd.append(['vg', 'annotate', '-p', '-x', os.path.basename(xg_file), '-a', '-'])
        cmd.append(['tee', os.path.basename(gam_annot_file)])
        cmd.append(['vg', 'view', '-aj', '-'])
        with open(gam_annot_json, 'w') as output_annot_json:
            options.drunner.call(job, cmd, work_dir = work_dir, outfile=output_annot_json)

        # turn the annotated gam json into truth positions, as separate command since
        # we're going to use a different docker container.  (Note, would be nice to
        # avoid writing the json to disk)        
        jq_cmd = ['jq', '-c', '-r', '[ .name, .refpos[0].name, .refpos[0].offset ] | @tsv',
                  os.path.basename(gam_annot_json)]

        # output truth positions
        true_pos_file = os.path.join(work_dir, 'true_{}.pos'.format(chunk_i))
        with open(true_pos_file, 'w') as out_true_pos:
            options.drunner.call(job, jq_cmd, work_dir = work_dir, outfile=out_true_pos)

        # get rid of that big json asap
        os.remove(gam_annot_json)

        # write to store. todo: there's probably no reason outside debugging to
        # keep both gams around.
        gam_chunk_id = write_to_store(job, options, gam_file)
        annot_gam_chunk_id = write_to_store(job, options, gam_annot_file)
        true_pos_chunk_id = write_to_store(job, options, true_pos_file)

        # return everythin as a tuple.
        return gam_chunk_id, annot_gam_chunk_id, true_pos_chunk_id
        

def run_merge_sim_chunks(job, options, sim_out_id_infos):
    """
    merge the sim output
    """
    assert len(sim_out_id_infos) > 0

    work_dir = job.fileStore.getLocalTempDir()

    if not options.gam:
        # merge up the reads files
        merged_reads_file = os.path.join('sim_reads')
        with open(merged_reads_file, 'a') as out_reads:
            for i, reads_file_id in enumerate(sim_out_id_infos):
                reads_file = os.path.join(work_dir, 'sim_reads_{}'.format(i))
                read_from_store(job, options, reads_file_id, reads_file)
                with open(reads_file) as rf:
                    shutil.copyfileobj(rf, out_reads)

        reads_id =  write_to_store(job, options, merged_reads_file)

        # checkpoint to the output store
        if not options.force_outstore:
            write_to_store(job, options, merged_reads_file, use_out_store = True)

        return reads_id
    
    else:
        # merge up the gam files
        merged_gam_file = os.path.join(work_dir, 'sim.gam')
        merged_annot_gam_file = os.path.join(work_dir, 'sim_annot.gam')
        merged_true_file = os.path.join(work_dir, 'true.pos.unsorted')
        
        with open(merged_gam_file, 'a') as out_gam, \
             open(merged_annot_gam_file, 'a') as out_annot_gam, \
             open(merged_true_file, 'a') as out_true:
            
            for i, sim_out_id_info in enumerate(sim_out_id_infos):
                gam_file = os.path.join(work_dir, 'sim_{}.gam'.format(i))
                read_from_store(job, options, sim_out_id_info[0], gam_file)
                with open(gam_file) as rf:
                    shutil.copyfileobj(rf, out_gam)
                    
                gam_annot_file = os.path.join(work_dir, 'sim_annot_{}.gam'.format(i))
                read_from_store(job, options, sim_out_id_info[1], gam_annot_file)
                with open(gam_annot_file) as rf:
                    shutil.copyfileobj(rf, out_annot_gam)

                true_file = os.path.join(work_dir, 'true_{}.pos'.format(i))
                read_from_store(job, options, sim_out_id_info[2], true_file)
                with open(true_file) as rf:
                    shutil.copyfileobj(rf, out_true)

        # sort the positions file
        sorted_true_file = os.path.join(work_dir, 'true.pos')
        sort_cmd = ['sort', os.path.basename(merged_true_file)]
        with open(sorted_true_file, 'w') as out_true:
            options.drunner.call(job, sort_cmd, work_dir = work_dir, outfile = out_true)

        
        merged_gam_id = write_to_store(job, options, merged_gam_file)
        merged_gam_annot_id = write_to_store(job, options, merged_annot_gam_file)
        true_id = write_to_store(job, options, sorted_true_file)

        # checkpoint to the output store
        if not options.force_outstore:
            write_to_store(job, options, merged_gam_file, use_out_store = True)
            write_to_store(job, options, merged_annot_gam_file, use_out_store = True)
            write_to_store(job, options, sorted_true_file, use_out_store = True)

        return merged_gam_id, merged_gam_annot_id, true_id
            
def sim_main(options):
    """
    Wrapper for vg sim. 
    """

    # make the docker runner
    options.drunner = ContainerRunner(
        container_tool_map = get_container_tool_map(options))

    require(all([i not in options.sim_opts for i in ['-x', '-n', '-a', '-s']]),
            ' sim-opts cannot contain -x, -n, -s or -a')
    require(options.sim_chunks > 0, '--sim_chunks must be >= 1')
    
    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'sim'

    # Throw error if something wrong with IOStore string
    IOStore.get(options.out_store)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            inputXGFileID = import_to_store(toil, options, options.xg_index)

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_sim, options,inputXGFileID,
                                     cores=options.misc_cores,
                                     memory=options.misc_mem,
                                     disk=options.misc_disk)
            
            # Run the job and store the returned list of output files to download
            toil.start(root_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
