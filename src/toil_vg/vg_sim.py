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
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_construct import run_unzip_fasta

logger = logging.getLogger(__name__)

def sim_subparser(parser):
    """
    Create a subparser for mapping.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("xg_indexes", nargs='+', type=make_url,
                        help="Path(s) to xg index(es) (separated by space)")
    parser.add_argument("num_reads", type=int,
                        help="Number of reads to simulate")
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--gam", action="store_true",
                        help="Output GAM file, annotated gam file and truth positions")
    parser.add_argument("--annotate_xg", type=make_url,
                        help="xg index used for gam annotation (if different from input indexes)")
    parser.add_argument("--sim_opts", type=str,
                        help="arguments for vg sim (wrapped in \"\"). Do not include -x, -n, -s, or -a")
    parser.add_argument("--sim_chunks", type=int, default=1,
                        help="split simulation into this many chunks, run in parallel when possible")
    parser.add_argument("--seed", type=int, default=None,
                        help="random seed")
    parser.add_argument("--fastq", type=make_url,
                        help="use error profile derived from given reads (ignores -l, -N, -f from sim_opts)")
                        
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def validate_sim_options(options):
    require(not options.sim_opts or all([i not in options.sim_opts for i in ['-x', '-n', '-a', '-s']]),
            ' sim-opts cannot contain -x, -n, -s or -a')
    require(options.sim_chunks > 0, '--sim_chunks must be >= 1')
    
def run_sim(job, context, num_reads, gam, seed, sim_chunks, xg_file_ids, xg_annot_file_id, fastq_id = None):
    """  
    run a bunch of simulation child jobs, merge up their output as a follow on
    """
    sim_out_id_infos = []

    # no seed specified, we choose one at random
    if seed is None:
        seed = random.randint(0, 2147483647)
        RealtimeLogger.info('No seed specifed, choosing random value = {}'.format(seed))

    # we can have more than one xg file if we've split our input graphs up
    # into haplotypes
    for xg_i, xg_file_id in enumerate(xg_file_ids):
        file_reads = num_reads / len(xg_file_ids)
        if xg_file_id == xg_file_ids[-1]:
            file_reads += num_reads % len(xg_file_ids)
        file_seed = seed + xg_i * sim_chunks
        
        # each element is either reads_chunk_id or (gam_chunk_id, annot_gam_chunk_id, true_pos_chunk_id)
        # if --gam not specified
        for chunk_i in range(sim_chunks):
            chunk_reads = file_reads / sim_chunks
            if chunk_i == sim_chunks - 1:
                chunk_reads += file_reads % sim_chunks
            sim_out_id_info = job.addChildJobFn(run_sim_chunk, context, gam, file_seed, xg_file_id,
                                                xg_annot_file_id,
                                                chunk_i, chunk_reads,
                                                fastq_id,
                                                cores=context.config.sim_cores, memory=context.config.sim_mem,
                                                disk=context.config.sim_disk).rv()
            sim_out_id_infos.append(sim_out_id_info)

    return job.addFollowOnJobFn(run_merge_sim_chunks, context, gam, sim_out_id_infos,
                                cores=context.config.sim_cores, memory=context.config.sim_mem,
                                disk=context.config.sim_disk).rv()

def run_sim_chunk(job, context, gam, seed, xg_file_id, xg_annot_file_id, chunk_i, num_reads, fastq_id):
    """
    simulate some reads (and optionally gam),
    return either reads_chunk_id or (gam_chunk_id, annot_gam_chunk_id, true_pos_chunk_id)
    if --gam specified
    """

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # read the xg file
    xg_file = os.path.join(work_dir, 'index.xg')
    job.fileStore.readGlobalFile(xg_file_id, xg_file)

    # read the fastq file
    if fastq_id:
        fastq_file = os.path.join(work_dir, 'error_template.fastq')
        job.fileStore.readGlobalFile(fastq_id, fastq_file)
    
    # and the annotation xg file
    if xg_annot_file_id:
        xg_annot_file = os.path.join(work_dir, 'annot_index.xg')
        xg_annot_file = job.fileStore.readGlobalFile(xg_annot_file_id, xg_annot_file)
    else:
        xg_annot_file = xg_file

    # run vg sim
    sim_cmd = ['vg', 'sim', '-x', os.path.basename(xg_file), '-n', num_reads] + context.config.sim_opts
    if seed is not None:
        sim_cmd += ['-s', seed + chunk_i]
    if fastq_id:
        sim_cmd += ['-F', os.path.basename(fastq_file)]

    if not gam:
        # output reads
        reads_file = os.path.join(work_dir, 'sim_reads_{}'.format(chunk_i))

        # run vg sim
        with open(reads_file, 'w') as output_reads:
            context.runner.call(job, sim_cmd, work_dir = work_dir, outfile=output_reads)

        # write to the store
        return context.write_intermediate_file(job, reads_file)
    else:
        # output gam
        gam_file = os.path.join(work_dir, 'sim_{}.gam'.format(chunk_i))
        gam_annot_file = os.path.join(work_dir, 'sim_{}_annot.gam'.format(chunk_i))
        gam_annot_json = os.path.join(work_dir, 'sim_{}_annot.json'.format(chunk_i))

        # run vg sim, write output gam, annotated gam, annotaged gam json
        # (from vg/scripts/map-sim)
        cmd = [sim_cmd + ['-a']]
        cmd.append(['tee', os.path.basename(gam_file)])
        cmd.append(['vg', 'annotate', '-p', '-x', os.path.basename(xg_annot_file), '-a', '-'])
        cmd.append(['tee', os.path.basename(gam_annot_file)])
        cmd.append(['vg', 'view', '-aj', '-'])
        with open(gam_annot_json, 'w') as output_annot_json:
            context.runner.call(job, cmd, work_dir = work_dir, outfile=output_annot_json)

        # turn the annotated gam json into truth positions, as separate command since
        # we're going to use a different docker container.  (Note, would be nice to
        # avoid writing the json to disk)        
        jq_cmd = ['jq', '-c', '-r', '[ .name, .refpos[0].name, .refpos[0].offset ] | @tsv',
                  os.path.basename(gam_annot_json)]

        # output truth positions
        true_pos_file = os.path.join(work_dir, 'true_{}.pos'.format(chunk_i))
        with open(true_pos_file, 'w') as out_true_pos:
            context.runner.call(job, jq_cmd, work_dir = work_dir, outfile=out_true_pos)

        # get rid of that big json asap
        os.remove(gam_annot_json)

        # write to store. todo: there's probably no reason outside debugging to
        # keep both gams around.
        gam_chunk_id = context.write_intermediate_file(job, gam_file)
        annot_gam_chunk_id = context.write_intermediate_file(job, gam_annot_file)
        true_pos_chunk_id = context.write_intermediate_file(job, true_pos_file)

        # return everythin as a tuple.
        return gam_chunk_id, annot_gam_chunk_id, true_pos_chunk_id
        

def run_merge_sim_chunks(job, context, gam, sim_out_id_infos):
    """
    merge the sim output
    """
    assert len(sim_out_id_infos) > 0

    work_dir = job.fileStore.getLocalTempDir()

    if not gam:
        # merge up the reads files
        merged_reads_file = os.path.join('sim_reads')
        with open(merged_reads_file, 'a') as out_reads:
            for i, reads_file_id in enumerate(sim_out_id_infos):
                reads_file = os.path.join(work_dir, 'sim_reads_{}'.format(i))
                job.fileStore.readGlobalFile(reads_file_id, reads_file)
                with open(reads_file) as rf:
                    shutil.copyfileobj(rf, out_reads)

        return context.write_output_file(job, merged_reads_file)
    
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
                job.fileStore.readGlobalFile(sim_out_id_info[0], gam_file)
                with open(gam_file) as rf:
                    shutil.copyfileobj(rf, out_gam)
                    
                gam_annot_file = os.path.join(work_dir, 'sim_annot_{}.gam'.format(i))
                job.fileStore.readGlobalFile(sim_out_id_info[1], gam_annot_file)
                with open(gam_annot_file) as rf:
                    shutil.copyfileobj(rf, out_annot_gam)

                true_file = os.path.join(work_dir, 'true_{}.pos'.format(i))
                job.fileStore.readGlobalFile(sim_out_id_info[2], true_file)
                with open(true_file) as rf:
                    shutil.copyfileobj(rf, out_true)

        # sort the positions file
        sorted_true_file = os.path.join(work_dir, 'true.pos')
        sort_cmd = ['sort', os.path.basename(merged_true_file)]
        with open(sorted_true_file, 'w') as out_true:
            context.runner.call(job, sort_cmd, work_dir = work_dir, outfile = out_true)

        
        merged_gam_id = context.write_output_file(job, merged_gam_file)
        merged_gam_annot_id = context.write_output_file(job, merged_annot_gam_file)
        true_id = context.write_output_file(job, sorted_true_file)

        return merged_gam_id, merged_gam_annot_id, true_id
            
def sim_main(context, options):
    """
    Wrapper for vg sim. 
    """

    validate_sim_options(options)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            inputXGFileIDs = []
            for xg_index in options.xg_indexes:
                inputXGFileIDs.append(toil.importFile(xg_index))
            if options.annotate_xg:
                inputAnnotXGFileID = toil.importFile(options.annotate_xg)
            else:
                inputAnnotXGFileID = None
            if options.fastq:
                inputFastqFileID = toil.importFile(options.fastq)
            else:
                inputFastqFileID = None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context)

            # Unzip the fastq
            if options.fastq and options.fastq.endswith('.gz'):
                inputFastqFileID = init_job.addChildJobFn(run_unzip_fasta, context, inputFastqFileID, 
                                                          os.path.basename(options.fastq)).rv()

            # Make a root job
            root_job = Job.wrapJobFn(run_sim, context, options.num_reads, options.gam,
                                     options.seed, options.sim_chunks,
                                     inputXGFileIDs,
                                     inputAnnotXGFileID,
                                     inputFastqFileID,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

            init_job.addFollowOn(root_job)            
            
            # Run the job and store the returned list of output files to download
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
