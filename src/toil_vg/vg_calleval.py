#!/usr/bin/env python2.7
"""
vg_calleval.py: Compare vcfs with vcfeval.  Option to make freebayes calls to use as baseline.  Can
run on vg_mapeval.py output. 
"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string, math
import urlparse
import getpass
import pdb
import gzip
import logging
import copy
from collections import Counter

from math import ceil
from subprocess import Popen, PIPE

try:
    import numpy as np
    from sklearn.metrics import roc_auc_score, average_precision_score, r2_score, roc_curve
    have_sklearn = True
except:
    have_sklearn = False

import tsv
import vcf

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.vg_call import chunked_call_parse_args, run_all_calling
from toil_vg.vg_vcfeval import vcfeval_parse_args, run_vcfeval, run_vcfeval_roc_plot
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_construct import run_unzip_fasta

logger = logging.getLogger(__name__)

def calleval_subparser(parser):
    """
    Create a subparser for calleval.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # Add the out_store
    # TODO: do this at a higher level?
    # Or roll into Context?
    parser.add_argument('out_store',
                        help='output store.  All output written here. Path specified using same syntax as toil jobStore')

    parser.add_argument("--chroms", nargs='+', required=True,
                        help="Name(s) of reference path in graph(s) (separated by space)."
                        " Must be same length/order as --gams")
    # todo: move to chunked_call_parse_args and share with toil-vg run
    parser.add_argument("--gams", nargs='+', required=True, type=make_url,
                        help="GAMs to call.  One per chromosome. Must be same length/order as --chroms")
    parser.add_argument("--sample_name", type=str, required=True,
                        help="sample name (ex NA12878)")
    parser.add_argument("--gam_index_cores", type=int,
                        help="number of threads used for gam indexing")


    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common call options shared with toil_vg pipeline
    chunked_call_parse_args(parser)
    
    # Add common vcfeval options shared with toil_vg pipeline
    vcfeval_parse_args(parser)

    # Add common calleval options shared with toil_vg pipeline
    calleval_parse_args(parser)

    # Add common docker options shared with toil_vg pipeline
    add_container_tool_parse_args(parser)
    
def calleval_parse_args(parser):
    """
    Add the calleval options to the given argparse parser.
    """
    parser.add_argument('--gam_names', nargs='+', required=True,
                        help='names of vg runs (corresponds to gams and xg_paths)')
    parser.add_argument('--xg_paths', nargs='+', required=True, type=make_url,
                        help='xg indexes for the different graphs')
    parser.add_argument('--freebayes', action='store_true',
                        help='run freebayes as a baseline')
    parser.add_argument('--bam_names', nargs='+',
                        help='names of bwa runs (corresponds to bams)')
    parser.add_argument('--bams', nargs='+', type=make_url,
                        help='bam inputs for freebayes')
    parser.add_argument('--call_and_genotype', action='store_true',
                        help='run both vg call and vg genotype')
    parser.add_argument('--filter_opts_gt',
                        help='override filter-opts for genotype only')
        
def validate_calleval_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(len(options.gam_names) == len(options.xg_paths) == len(options.gams),
            '--gam_names, --xg_paths, --gams must all contain same number of elements')
    if options.freebayes:
        require(options.bams, '--bams must be given for use with freebayes')
    if options.bams or options.bam_names:
        require(options.bams and options.bam_names and len(options.bams) == len(options.bam_names),
                '--bams and --bam_names must be same length')
    # todo: generalize.  
    require(options.chroms and len(options.chroms) == 1,
            'one sequence must be specified with --chroms')
    require(options.vcfeval_baseline, '--vcfeval_baseline required')
    require(options.vcfeval_fasta, '--vcfeval_fasta required')

def run_bam_index(job, context, bam_file_id, bam_name):
    """
    sort and index a bam.  return the sorted bam and its idx
    """
    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download the input
    bam_path = os.path.join(work_dir, bam_name + '.bam')
    job.fileStore.readGlobalFile(bam_file_id, bam_path)

    sort_bam_path = os.path.join(work_dir, 'sort.bam')
    sort_cmd = ['samtools', 'sort', os.path.basename(bam_path), '-o',
                os.path.basename(sort_bam_path), '-O', 'BAM', '-@', max(0, job.cores - 1)]
    context.runner.call(job, sort_cmd, work_dir=work_dir)
    bam_index_cmd = ['samtools', 'index', os.path.basename(sort_bam_path)]
    context.runner.call(job, bam_index_cmd, work_dir=work_dir)

    out_bam_id = context.write_intermediate_file(job, sort_bam_path)
    out_idx_id = context.write_intermediate_file(job, sort_bam_path + '.bai')
    return out_bam_id, out_idx_id
    
    
def run_freebayes(job, context, fasta_file_id, bam_file_id, bam_idx_id,
                  sample_name, region, offset, out_name,
                  freebayes_opts = ['--genotype-qualities']):
    """
    run freebayes to make a vcf
    """

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download the input
    fasta_path = os.path.join(work_dir, 'ref.fa')
    bam_path = os.path.join(work_dir, 'alignment.bam')
    bam_idx_path = bam_path + '.bai'
    job.fileStore.readGlobalFile(fasta_file_id, fasta_path)
    job.fileStore.readGlobalFile(bam_file_id, bam_path)
    job.fileStore.readGlobalFile(bam_idx_id, bam_idx_path)

    # run freebayes
    fb_cmd = ['freebayes', '-f', os.path.basename(fasta_path), os.path.basename(bam_path)]
    if freebayes_opts:
        fb_cmd += freebayes_opts

    if region:
        fb_cmd += ['-r', region]

    vcf_path = os.path.join(work_dir, '{}-raw.vcf'.format(out_name))
    timer = TimeTracker('freebayes')
    with open(vcf_path, 'w') as out_vcf:
        context.runner.call(job, fb_cmd, work_dir=work_dir, outfile=out_vcf)
    timer.stop()

    context.write_intermediate_file(job, vcf_path)

    vcf_fix_path = os.path.join(work_dir, '{}.vcf'.format(out_name))
    
    # apply offset and sample name
    vcf_reader = vcf.Reader(open(vcf_path))
    vcf_writer = vcf.Writer(open(vcf_fix_path, 'w'), vcf_reader)
    for record in vcf_reader:
        if offset:
            record.POS += int(offset)
        if sample_name:
            pass
        vcf_writer.write_record(record)
    vcf_writer.flush()
    vcf_writer.close()

    context.runner.call(job, ['bgzip', os.path.basename(vcf_fix_path)], work_dir = work_dir)
    context.runner.call(job, ['tabix', '-p', 'vcf', os.path.basename(vcf_fix_path) + '.gz'], work_dir = work_dir)

    return (context.write_output_file(job, vcf_fix_path + '.gz'),
            context.write_output_file(job, vcf_fix_path + '.gz.tbi'),
            timer)

def run_calleval_results(job, context, names, vcf_tbi_pairs, eval_results, timing_results):
    """ output the calleval results"""

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    has_clipped_results = all([eval_result[1] is not None for eval_result in eval_results])
    result_idx = 1 if has_clipped_results else 0

    # make a simple tsv
    stats_path = os.path.join(work_dir, 'calleval_stats.tsv')
    with open(stats_path, 'w') as stats_file:
        for name, eval_result, in zip(names, eval_results):
            stats_file.write('{}\t{}\n'.format(name, eval_result[result_idx][0]))

    # make some roc plots
    roc_plot_ids = []
    for i, roc_type in zip(range(3,6), ['snp', 'non_snp', 'weighted']):
        roc_table_ids = [eval_result[0][i] for eval_result in eval_results]
        roc_title = roc_type if not has_clipped_results else roc_type + '-unclipped'
        roc_plot_ids.append(job.addChildJobFn(run_vcfeval_roc_plot, context, roc_table_ids, names=names,
                                              title=roc_title).rv())
        
        if has_clipped_results:
            roc_table_clip_ids = [eval_result[1][i] for eval_result in eval_results]
            roc_plot_ids.append(job.addChildJobFn(run_vcfeval_roc_plot, context, roc_table_clip_ids, names=names,
                                                  title=roc_type).rv())

    # write some times
    times_path = os.path.join(work_dir, 'call_times.tsv')
    
    # organize our expected labels a bit
    calling_labels = ['call', 'genotype', 'freebayes']
    augmenting_labels = ['call-filter-augment', 'call-augment', 'call-filter']
    all_labels = set()
    for timer in timing_results:
        for label in timer.names():
            all_labels.add(label)
    other_labels = all_labels - set(calling_labels + augmenting_labels)
    with open(times_path, 'w') as times_file:
        # write the header
        times_file.write('method\tcall\taugment\ttotal-other')
        for other_label in other_labels:
            times_file.write('\t{}'.format(other_label.replace('call-','')))
        times_file.write('\n')
        # write the body
        for name, timer in zip(names, timing_results):
            times_file.write('{}\t{}'.format(name, timer.total(calling_labels)))
            times_file.write('\t{}'.format(timer.total(augmenting_labels)))
            times_file.write('\t{}'.format(timer.total(other_labels)))
            for other_label in other_labels:
                times_file.write('\t{}'.format(timer.total([other_label])))
            times_file.write('\n')
                        
    return [context.write_output_file(job, stats_path),
            context.write_output_file(job, times_path)] + roc_plot_ids
                             
        
def run_calleval(job, context, xg_ids, gam_ids, bam_ids, bam_idx_ids, gam_names, bam_names,
                 vcfeval_baseline_id, vcfeval_baseline_tbi_id, fasta_id, bed_id,
                 call, genotype, sample_name, chrom, vcf_offset, vcfeval_score_field,
                 filter_opts_gt):
    """ top-level call-eval function.  runs the caller and genotype on every gam,
    and freebayes on every bam.  the resulting vcfs are put through vcfeval
    and the accuracies are tabulated in the output
    """
    vcf_tbi_id_pairs = [] 
    names = []
    eval_results = []
    timing_results = []

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)
    
    if bam_ids:
        for bam_id, bam_idx_id, bam_name in zip(bam_ids, bam_idx_ids, bam_names):
            if not bam_idx_id:
                bam_index_job = child_job.addChildJobFn(run_bam_index, context, bam_id, bam_name,
                                                        cores=context.config.calling_cores,
                                                        memory=context.config.calling_mem,
                                                        disk=context.config.calling_disk)
                sorted_bam_id = bam_index_job.rv(0)
                sorted_bam_idx_id = bam_index_job.rv(1)
            else:
                bam_index_job = Job()
                child_job.addChild(bam_index_job)
                sorted_bam_id = bam_id
                sorted_bam_idx_id = bam_idx_id                

            fb_out_name = '{}-fb'.format(bam_name)
            fb_job = bam_index_job.addFollowOnJobFn(run_freebayes, context, fasta_id,
                                                    sorted_bam_id, sorted_bam_idx_id, sample_name,
                                                    chrom, vcf_offset,
                                                    out_name = fb_out_name,
                                                    cores=context.config.calling_cores,
                                                    memory=context.config.calling_mem,
                                                    disk=context.config.calling_disk)
            fb_vcf_tbi_id_pair = (fb_job.rv(0), fb_job.rv(1))
            timing_results.append(fb_job.rv(2))

            if bed_id:
                eval_clip_result = fb_job.addFollowOnJobFn(run_vcfeval, context, sample_name, fb_vcf_tbi_id_pair,
                                                           vcfeval_baseline_id, vcfeval_baseline_tbi_id, 'ref.fasta',
                                                           fasta_id, bed_id, out_name=fb_out_name,
                                                           score_field='GQ').rv()
            else:
                eval_clip_result = None
                
            eval_result = fb_job.addFollowOnJobFn(run_vcfeval, context, sample_name, fb_vcf_tbi_id_pair,
                                                  vcfeval_baseline_id, vcfeval_baseline_tbi_id, 'ref.fasta',
                                                  fasta_id, None,
                                                  out_name=fb_out_name if not bed_id else fb_out_name + '-unclipped',
                                                  score_field='GQ').rv()
            
            vcf_tbi_id_pairs.append(fb_vcf_tbi_id_pair)            
            names.append(fb_out_name)
            eval_results.append((eval_result, eval_clip_result))

    # optional override to filter-opts when running genotype
    # this is allows us to run different filter-opts for call and genotype
    # (something we don't really need an interface for besides in calleval)
    if filter_opts_gt:
        gt_context = copy.deepcopy(context)
        gt_context.config.filter_opts = filter_opts_gt.split()
    else:
        gt_context = context

    if gam_ids:
        for gam_id, gam_name, xg_id in zip(gam_ids, gam_names, xg_ids):
            for gt in [False, True]:
                if (call and not gt) or (genotype and gt):
                    out_name = '{}{}'.format(gam_name, '-gt' if gt else '-call')
                    call_job = child_job.addChildJobFn(run_all_calling, gt_context if gt else context,
                                                       xg_id, [gam_id], [chrom], [vcf_offset],
                                                       sample_name, genotype=gt,
                                                       out_name=out_name,
                                                       cores=context.config.misc_cores,
                                                       memory=context.config.misc_mem,
                                                       disk=context.config.misc_disk)
                    vcf_tbi_id_pair = (call_job.rv(0), call_job.rv(1))
                    timing_results.append(call_job.rv(2))

                    if not vcfeval_score_field:
                        score_field = 'GQ' if gt else 'QUAL'
                    else:
                        score_field = vcfeval_score_field

                    if bed_id:
                        eval_clip_result = call_job.addFollowOnJobFn(run_vcfeval, context, sample_name, vcf_tbi_id_pair,
                                                                     vcfeval_baseline_id, vcfeval_baseline_tbi_id, 'ref.fasta',
                                                                     fasta_id, bed_id, out_name=out_name,
                                                                     score_field=score_field).rv()
                    else:
                        eval_clip_result = None
                        
                    eval_result = call_job.addFollowOnJobFn(run_vcfeval, context, sample_name, vcf_tbi_id_pair,
                                                            vcfeval_baseline_id, vcfeval_baseline_tbi_id, 'ref.fasta',
                                                            fasta_id, None,
                                                            out_name=out_name if not bed_id else out_name + '-unclipped',
                                                            score_field=score_field).rv()
                        
                    names.append(out_name)            
                    vcf_tbi_id_pairs.append(vcf_tbi_id_pair)
                    eval_results.append((eval_result, eval_clip_result))

    calleval_results = child_job.addFollowOnJobFn(run_calleval_results, context, names,
                                                  vcf_tbi_id_pairs, eval_results, timing_results,
                                                  cores=context.config.misc_cores,
                                                  memory=context.config.misc_mem,
                                                  disk=context.config.misc_disk).rv()

    return calleval_results, names, vcf_tbi_id_pairs, eval_results

def calleval_main(context, options):
    """ entrypoint for calling """

    validate_calleval_options(options)
            
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()

            # Upload local files to the job store            
            inputXGFileIDs = []
            xgToID = {}
            if options.xg_paths:
                for xg_path in options.xg_paths:
                    # we allow same files to be passed many times, but just import them once                    
                    if xg_path not in xgToID:
                        xgToID[xg_path] = toil.importFile(xg_path)
                    inputXGFileIDs.append(xgToID[xg_path])
            inputGamFileIDs = []
            gamToID = {}
            if options.gams:
                for gam in options.gams:
                    if gam not in gamToID:
                        gamToID[gam] = toil.importFile(gam)
                    inputGamFileIDs.append(gamToID[gam])
                        
            inputBamFileIDs = []
            inputBamIdxIds = []
            if options.bams:
                for bam in options.bams:
                    inputBamFileIDs.append(toil.importFile(bam))
                    try:
                        bamIdxId = toil.importFile(bam + '.bai')
                    except:
                        bamIdxId = None
                    inputBamIdxIds.append(bamIdxId)

            vcfeval_baseline_id = toil.importFile(options.vcfeval_baseline)
            vcfeval_baseline_tbi_id = toil.importFile(options.vcfeval_baseline + '.tbi')
            fasta_id = toil.importFile(options.vcfeval_fasta)
            bed_id = toil.importFile(options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)

            if options.vcfeval_fasta.endswith('.gz'):
                # unzip the fasta
                fasta_id = init_job.addChildJobFn(run_unzip_fasta, context, fasta_id,
                                                  os.path.basename(options.vcfeval_fasta)).rv()

            # Make a root job
            do_call = options.call_and_genotype or not options.genotype
            do_genotype = options.call_and_genotype or options.genotype
            root_job = Job.wrapJobFn(run_calleval, context, inputXGFileIDs, inputGamFileIDs, inputBamFileIDs,
                                     inputBamIdxIds,
                                     options.gam_names, options.bam_names, 
                                     vcfeval_baseline_id, vcfeval_baseline_tbi_id, fasta_id, bed_id,
                                     do_call,
                                     do_genotype,
                                     options.sample_name,
                                     options.chroms[0], options.vcf_offsets[0] if options.vcf_offsets else 0,
                                     options.vcfeval_score_field,
                                     options.filter_opts_gt,
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
    
    
    
