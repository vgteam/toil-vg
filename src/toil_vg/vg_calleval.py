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
from toil_vg.vg_call import chunked_call_parse_args, run_all_calling, run_merge_vcf
from toil_vg.vg_vcfeval import vcfeval_parse_args, run_extract_sample_truth_vcf, run_vcfeval, run_vcfeval_roc_plot, run_happy, run_sv_eval
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_construct import run_unzip_fasta
from toil_vg.vg_surject import run_surjecting

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
                        help="Name(s) of reference path in graph(s) (separated by space).")
    # todo: move to chunked_call_parse_args and share with toil-vg run
    parser.add_argument("--gams", nargs='+', type=make_url,
                        help="GAMs to call.  Each GAM treated as separate input (and must contain all chroms)")
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
    parser.add_argument('--gam_names', nargs='+',
                        help='names of vg runs (corresponds to gams and xg_paths)')
    parser.add_argument('--xg_paths', nargs='+', type=make_url,
                        help='xg indexes for the different graphs')
    parser.add_argument('--freebayes', action='store_true',
                        help='run freebayes as a baseline')
    parser.add_argument("--caller_fasta", type=make_url,
                        help="Use this FASTA instead of the vcfeval fasta for Freebayes. Maybe be gzipped")
    parser.add_argument('--platypus', action='store_true',
                        help='run platypus as a baseline')
    parser.add_argument('--bam_names', nargs='+', default=[],
                        help='names of bwa runs (corresponds to bams)')
    parser.add_argument('--bams', nargs='+', type=make_url, default=[],
                        help='bam inputs for freebayes')
    parser.add_argument('--filter_opts_gt',
                        help='override filter-opts for genotype only')
    parser.add_argument('--clip_only', action='store_true',
                        help='only compute accuracy clipped to --vcfeval_bed_regions')
    parser.add_argument('--plot_sets', nargs='+', default=[],
                        help='comma-separated lists of condition-tagged BAM/GAM names (primary-mp-pe, etc.) with colon-separated title prefixes')
    parser.add_argument('--call', action='store_true',
                        help='run vg call (even if --genotype used)')    
    parser.add_argument("--surject", action="store_true",
                        help="surject GAMs to BAMs, adding the latter to the comparison")
    parser.add_argument("--interleaved", action="store_true", default=False,
                        help="assume GAM files are interleaved when surjecting")
        
def validate_calleval_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    if options.gams or options.gam_names or options.xg_paths:        
        require(options.gams and options.gam_names and options.xg_paths and
                len(options.gam_names) == len(options.xg_paths) == len(options.gams),
                '--gam_names, --xg_paths, --gams must all contain same number of elements')
        require(options.call or options.genotype or options.surject,
                '--call and/or --genotype and/or --surject required with --gams')
    if options.freebayes or options.platypus:
        require(options.bams or options.surject, '--bams and/or --surject needed with --freebayes or --platypus')
    if options.bams or options.bam_names:
        require(options.bams and options.bam_names and len(options.bams) == len(options.bam_names),
                '--bams and --bam_names must be same length')
    require(options.vcfeval_baseline, '--vcfeval_baseline required')
    require(options.vcfeval_fasta, '--vcfeval_fasta required')
    if options.call or options.genotype or options.surject:
        require(options.gams, '--gams must be given with --call, --genotype, and --surject')
    require(options.vcfeval_bed_regions is not None or not options.clip_only,
            '--vcfeval_bed_regions must be given with --clip_only')
    if options.surject or options.bams:
        require(options.freebayes or options.platypus,
                '--freebayes and/or --platypus must be used with --bams and --surject')

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
    
def run_all_bam_caller(job, context, fasta_file_id, bam_file_id, bam_idx_id,
                       sample_name, chroms, offsets, out_name, bam_caller,
                       bam_caller_opts = []):
    """
    run freebayes or platypus on a set of chromosomal regions.  this is done by sending each region to a 
    child job and farming off the entire input to each (ie not splitting the input)
    """
    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)

    fb_vcf_ids = []
    fb_tbi_ids = []
    fb_timers = []
    assert chroms
    if not offsets:
        offsets = [None] * len(chroms)
    for chrom, offset in zip(chroms, offsets):
        fb_job = child_job.addChildJobFn(run_bam_caller, context, fasta_file_id, bam_file_id, bam_idx_id,
                                         sample_name, chrom, offset, out_name, bam_caller, bam_caller_opts,
                                         cores=context.config.calling_cores,
                                         memory=context.config.calling_mem,
                                         disk=context.config.calling_disk)
        fb_vcf_ids.append(fb_job.rv(0))
        fb_tbi_ids.append(fb_job.rv(1))
        fb_timers.append([fb_job.rv(2)])

    merge_vcf_job = child_job.addFollowOnJobFn(run_merge_vcf, context, out_name, zip(fb_vcf_ids, fb_tbi_ids), fb_timers)
    return merge_vcf_job.rv()
    
def run_bam_caller(job, context, fasta_file_id, bam_file_id, bam_idx_id,
                   sample_name, chrom, offset, out_name, bam_caller, bam_caller_opts):
    """
    run freebayes or platypus to make a vcf
    """

    assert bam_caller in ['freebayes', 'platypus']
    
    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download the input
    fasta_path = os.path.join(work_dir, 'ref.fa')
    bam_path = os.path.join(work_dir, 'alignment.bam')
    bam_idx_path = bam_path + '.bai'
    job.fileStore.readGlobalFile(fasta_file_id, fasta_path)
    job.fileStore.readGlobalFile(bam_file_id, bam_path)
    job.fileStore.readGlobalFile(bam_idx_id, bam_idx_path)

    # output
    vcf_path = os.path.join(work_dir, '{}-raw.vcf'.format(out_name))
    
    if bam_caller == 'freebayes':
        fb_cmd = ['freebayes', '-f', os.path.basename(fasta_path), os.path.basename(bam_path)]
        if bam_caller_opts:
            fb_cmd += bam_caller_opts
        if chrom:
            fb_cmd += ['-r', chrom]
        timer = TimeTracker('freebayes')
        with open(vcf_path, 'w') as out_vcf:
            context.runner.call(job, fb_cmd, work_dir=work_dir, outfile=out_vcf)
        timer.stop()
            
    elif bam_caller == 'platypus':
        context.runner.call(job, ['samtools', 'faidx', 'ref.fa'], work_dir=work_dir)
        plat_cmd = ['Platypus.py', 'callVariants', '--refFile', os.path.basename(fasta_path),
                    '--bamFiles', os.path.basename(bam_path), '-o', os.path.basename(vcf_path)]
        if bam_caller_opts:
            plat_cmd += bam_caller_opts
        if chrom:
            plat_cmd += ['--regions', chrom]
        timer = TimeTracker('Platypus')
        context.runner.call(job, plat_cmd, work_dir=work_dir)
        timer.stop()
    else:
        assert False

    context.write_intermediate_file(job, vcf_path)

    vcf_fix_path = os.path.join(work_dir, '{}.vcf'.format(out_name))
    
    # apply offset and sample name
    vcf_reader = vcf.Reader(open(vcf_path))
    if sample_name:
        # Freebayes always outputs "unknown" for the sample if the BAM
        # doesn't specify a name, and can't be convinced to do otherwise:
        # https://github.com/ekg/freebayes/issues/471
        # So we hack the VCFReader to think the sample names are what we want them to be
        assert(len(vcf_reader.samples) == 1)
        # TODO: asserting that sample 0 is named unknown can fail. Why?
        RealtimeLogger.info('Correcting Freebayes samples {} to [{}]'.format(vcf_reader.samples, sample_name))
        vcf_reader.samples = [sample_name]
        # Rebuild the secret sample index
        vcf_reader._sample_indexes = {sample_name: 0}
    vcf_writer = vcf.Writer(open(vcf_fix_path, 'w'), vcf_reader)
    
    # We insist that Freebayes produce at least one record, or we will abort and dump files.
    have_records = False
    for record in vcf_reader:
        have_records = True
        if offset:
            record.POS += int(offset)
        vcf_writer.write_record(record)
    vcf_writer.flush()
    vcf_writer.close()
    
    if not have_records:
        # This is a problem. Stop here for debugging.
        RealtimeLogger.error('{} produced no calls. Dumping files...'.format(bam_caller))
        for dump_path in [fasta_path, bam_path, bam_idx_path]:
            if dump_path and os.path.isfile(dump_path):
                context.write_output_file(job, dump_path)
        raise RuntimeError('{} produced no calls'.format(bam_caller))

    context.runner.call(job, ['bgzip', os.path.basename(vcf_fix_path)], work_dir = work_dir)
    context.runner.call(job, ['tabix', '-p', 'vcf', os.path.basename(vcf_fix_path) + '.gz'], work_dir = work_dir)

    return (context.write_output_file(job, vcf_fix_path + '.gz'),
            context.write_output_file(job, vcf_fix_path + '.gz.tbi'),
            timer)

def run_calleval_plots(job, context, names, eval_results_dict, plot_sets):
    """
    
    Make and output calleval ROC plots.
    
    Takes a "names" list of all conditions. Condition names in the list (or in
    plot_sets) do not include an "-unclipped" tag; both clipped and unclipped
    plots are made if the data is available.
    
    Takes a nested dict by condition name, then clipping status ("clipped",
    "unclipped"), and then variant type ("snp", "non_snp", "weighted").
    Eventual entries are to ROC data file ids (.tsv.gz).
    
    plot_sets is a data structure of collections of conditions to plot against
    each other, as produced by parse_plot_sets.
    
    Returns a list of created plot file IDs.
    
    """
    
    # make some roc plots
    roc_plot_ids = []
    for roc_type in ['snp', 'non_snp', 'weighted']:
        # For each type of ROC to plot
    
        for mode in ['unclipped', 'clipped']:
            # For each clipping mode

            # What kind of ROC plot is it?
            # It should be this unless there is a clipped mode for this ROC and we are the unclipped one.
            roc_kind = roc_type
            if mode == 'unclipped' and True in [eval_results_dict.has_key('clipped') for eval_result in eval_results_dict.itervalues()]:
                # We are the unclipped mode and there will be a clipped mode
                roc_kind += '-unclipped'
            
            # Get all the eval results for this mode, by condition name
            mode_results = {name: eval_result.get(mode) for name, eval_result in eval_results_dict.iteritems()}
            
            if None in mode_results.itervalues():
                # We can't do this mode since it wasn't run
                continue
                
            # Extract out all the stats file IDs for this ROC type, by
            # condition name. This dict is implicitly for the clipped or
            # unclipped mode, depending on which mode we are doing. The
            # condition name keys don't have clipping tags.
            roc_table_ids = {name: result.get(roc_type) for name, result in mode_results.iteritems()}
            
            for subset_number, plot_set in enumerate(plot_sets):
                # For each plot set specifier
                
                # Unpack plot_set
                plot_title, plot_conditions = plot_set
                
                if plot_conditions is None:
                    # Use all the conditions instead
                    plot_conditions = names
                    
                # Add the ROC plot kind (snp, clipped, etc.) to the plot title
                if plot_title is None:
                    plot_title = roc_kind
                else:
                    plot_title = '{} ({})'.format(plot_title, roc_kind)
                
                for name in plot_conditions:
                    if name not in names:
                        # Complain if any of the names in the subset aren't conditions that actually ran
                        message = 'Condition {} not found in list of available conditions {}'.format(name, names)
                        RealtimeLogger.error(message)
                        raise RuntimeError(message)
                    if not roc_table_ids.has_key(name):
                        # Complain if any of the names in the subset lacks ROC data
                        message = 'Condition {} has no ROC data; data only available for {}'.format(name, roc_table_ids.keys())
                        RealtimeLogger.error(message)
                        raise RuntimeError(message)

                # Make a list of roc table IDs for the subset, in plot_conditions order
                subset_ids = [roc_table_ids[name] for name in plot_conditions]
                
                # Make the plot
                roc_plot_ids.append(job.addChildJobFn(run_vcfeval_roc_plot, context, subset_ids, names=plot_conditions,
                                                      kind=roc_kind, number=subset_number, title=plot_title).rv())
                                                      
    return roc_plot_ids
    

def run_calleval_results(job, context, names, vcf_tbi_pairs, eval_results_dict, happy_results_dict, sveval_results_dict,
                         timing_results, plot_sets):
    """
    
    output the calleval results
    
    Requires that, if any result in eval_results has clipped results, all
    results have clipped results, and similarly for unclipped results.

    plot_sets is a data structure of collections of conditions to plot against
    each other, as produced by parse_plot_sets.
    
    """

    RealtimeLogger.info('Handling results for conditions {} in sets {}'.format(names, plot_sets))

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # make a simple tsv
    stats_path = os.path.join(work_dir, 'calleval_stats.tsv')
    with open(stats_path, 'w') as stats_file:
        for name in names:
            eval_result = eval_results_dict[name]
            happy_result = happy_results_dict[name]
            # Find the best result (clipped if present, unclipped if not).
            best_result = eval_result.get("clipped", eval_result.get("unclipped", None))

            # Same for the happy results
            best_happy_result = happy_result.get("clipped", happy_result.get("unclipped", None))
            if best_happy_result is not None:
                happy_snp_f1 = best_happy_result['parsed_summary']['SNP']['METRIC.F1_Score']
                happy_indel_f1 =  best_happy_result['parsed_summary']['INDEL']['METRIC.F1_Score']
            else:
                happy_snp_f1, happy_indel_f1 = -1, -1

            # Same for the sveval results
            sveval_result = sveval_results_dict[name]
            best_sveval_result = sveval_result.get("clipped", sveval_result.get("unclipped", None))
            sveval_f1 = best_sveval_result['F1'] if best_sveval_result is not None else -1
                
            # Output the F1 scores
            # TODO: We amy want sveval_f1 in here, but the tests prohibit tis presence
            stats_file.write('{}\t{}\t{}\t{}\n'.format(name, best_result['f1'], happy_snp_f1, happy_indel_f1))

    # Make the roc plots
    roc_plot_job = job.addChildJobFn(run_calleval_plots, context, names, eval_results_dict, plot_sets)
    roc_plot_ids = roc_plot_job.rv()

    # write some times
    times_path = os.path.join(work_dir, 'call_times.tsv')
    
    # organize our expected labels a bit
    calling_labels = ['call', 'genotype', 'freebayes', 'platypus']
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
                        
    return roc_plot_job.addFollowOnJobFn(run_concat_lists,
                                         [context.write_output_file(job, stats_path), context.write_output_file(job, times_path)],
                                         roc_plot_ids).rv()
            
def run_vcf_subset(job, context, vcf_file_id, tbi_file_id, regions):
    # use bcftools to clip vcf to regions
    work_dir = job.fileStore.getLocalTempDir()

    # download the input
    vcf_path = os.path.join(work_dir, 'input.vcf.gz')
    out_vcf_path = os.path.join(work_dir, 'output.vcf.gz')
    job.fileStore.readGlobalFile(vcf_file_id, vcf_path)
    job.fileStore.readGlobalFile(tbi_file_id, vcf_path + '.tbi')

    cmd = ['bcftools', 'view', os.path.basename(vcf_path), '--regions', ','.join(regions),
           '--output-type', 'z', '--output-file', os.path.basename(out_vcf_path)]
    context.runner.call(job, cmd, work_dir=work_dir)
    context.runner.call(job, ['tabix', '--preset', 'vcf', os.path.basename(out_vcf_path)], work_dir=work_dir)
    return (context.write_intermediate_file(job, out_vcf_path),
            context.write_intermediate_file(job, out_vcf_path + '.tbi'))                             
        
def run_calleval(job, context, xg_ids, gam_ids, gam_idx_ids, bam_ids, bam_idx_ids, gam_names, bam_names,
                 vcfeval_baseline_id, vcfeval_baseline_tbi_id, caller_fasta_id, vcfeval_fasta_id,
                 bed_id, clip_only, call, genotype, sample_name, chroms, vcf_offsets,
                 vcfeval_score_field, plot_sets, filter_opts_gt, surject, interleaved,
                 freebayes, platypus, happy, sveval, recall):
    """
    top-level call-eval function. Runs the caller and genotype on every
    gam, and freebayes on every bam. The resulting vcfs are put through
    vcfeval and the accuracies are tabulated in the output.
    
    Returns the output of run_calleval results, a list of condition names, a
    list of corresponding called VCF.gz and index ID pairs, and dicts of
    vcfeval and happy result dicts, by condition name and clipped/unclipped
    status.

    plot_sets is a data structure of collections of conditions to plot against
    each other, as produced by parse_plot_sets.
    
    """
    
    # We store the name of each condition we run
    names = []
    # And we build up these result lists in sync with the name list
    vcf_tbi_id_pairs = []
    timing_results = []
    
    # Here we accumulate vcf_eval comparison results in a dict by condition name, then clipping status ("clipped", "unclipped").
    # Each contained dict is the output dict from run_vcfeval
    eval_results = collections.defaultdict(dict)
    # And here we store similarly the output dicts from run_happy
    happy_results = collections.defaultdict(dict)
    # And here we store similarly the output dicts from run_sveval
    sveval_results = collections.defaultdict(dict)

    # Some prep work (surjection and truth extraction) will happen under this head job
    head_job = Job()
    job.addChild(head_job)

    # Most of our work will run under this child job
    child_job = Job()
    head_job.addFollowOn(child_job)
    
    
    # We always extract a single-sample VCF from the truth, to save time
    # picking through all its samples multiple times over later. This should
    # also save memory. TODO: should we define a separate disk/memory requirement set?
    sample_extract_job = head_job.addChildJobFn(run_extract_sample_truth_vcf, context, sample_name, vcfeval_baseline_id,
                                                vcfeval_baseline_tbi_id, cores=context.config.vcfeval_cores,
                                                memory=context.config.vcfeval_mem,
                                                disk=context.config.vcfeval_disk)
    truth_vcf_id = sample_extract_job.rv(0)
    truth_vcf_tbi_id = sample_extract_job.rv(1)

    if not gam_idx_ids:
        gam_idx_ids = [None] * len(gam_ids)
    assert len(gam_idx_ids) == len(gam_ids)
    
    if surject:
        # optionally surject all the gams into bams
        for xg_id, gam_name, gam_id in zip(xg_ids, gam_names, gam_ids):
            surject_job = head_job.addChildJobFn(run_surjecting, context, gam_id, gam_name + '-surject',
                                                 interleaved, xg_id, chroms, cores=context.config.misc_cores,
                                                 memory=context.config.misc_mem, disk=context.config.misc_disk)
            bam_ids.append(surject_job.rv())
            bam_idx_ids.append(None)
            bam_names.append(gam_name + '-surject')
    
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

            bam_caller_infos = []
            if freebayes:
                bam_caller_infos.append(('freebayes', ['--genotype-qualities'], '-fb'))
            if platypus:
                bam_caller_infos.append(('platypus', ['--mergeClusteredVariants=1'], '-plat'))
                
            for bam_caller, bam_caller_opts, bam_caller_tag in bam_caller_infos:

                bam_caller_out_name = '{}{}'.format(bam_name, bam_caller_tag)
                bam_caller_job = bam_index_job.addFollowOnJobFn(run_all_bam_caller, context, caller_fasta_id,
                                                                sorted_bam_id, sorted_bam_idx_id, sample_name,
                                                                chroms, vcf_offsets,
                                                                out_name = bam_caller_out_name,
                                                                bam_caller = bam_caller,
                                                                bam_caller_opts = bam_caller_opts,
                                                                cores=context.config.misc_cores,
                                                                memory=context.config.misc_mem,
                                                                disk=context.config.misc_disk)
                bam_caller_vcf_tbi_id_pair = (bam_caller_job.rv(0), bam_caller_job.rv(1))
                timing_result = bam_caller_job.rv(2)

                if bed_id:

                    eval_results[bam_caller_out_name]["clipped"] = \
                        bam_caller_job.addFollowOnJobFn(run_vcfeval, context, sample_name, bam_caller_vcf_tbi_id_pair,
                                                        truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                        vcfeval_fasta_id, bed_id, out_name=bam_caller_out_name,
                                                        score_field='GQ', cores=context.config.vcfeval_cores,
                                                        memory=context.config.vcfeval_mem,
                                                        disk=context.config.vcfeval_disk).rv()
                    if happy:
                        happy_results[bam_caller_out_name]["clipped"] = \
                        bam_caller_job.addFollowOnJobFn(run_happy, context, sample_name, bam_caller_vcf_tbi_id_pair,
                                                        truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                        vcfeval_fasta_id, bed_id, out_name=bam_caller_out_name,
                                                        cores=context.config.vcfeval_cores,
                                                        memory=context.config.vcfeval_mem,
                                                        disk=context.config.vcfeval_disk).rv()

                    if sveval:
                        sveval_results[bam_caller_out_name]["clipped"] = \
                        bam_caller_job.addFollowOnJobFn(run_sv_eval, context, sample_name, bam_caller_vcf_tbi_id_pair,
                                                        truth_vcf_id, truth_vcf_tbi_id, min_sv_len=25,
                                                        sv_overlap=0.5, sv_region_overlap=1,
                                                        bed_id=bed_id,
                                                        out_name=bam_caller_out_name,
                                                        cores=context.config.vcfeval_cores,
                                                        memory=context.config.vcfeval_mem,
                                                        disk=context.config.vcfeval_disk).rv()                    

                if not clip_only:
                    # Also do unclipped

                    eval_results[bam_caller_out_name]["unclipped"] = \
                        bam_caller_job.addFollowOnJobFn(run_vcfeval, context, sample_name, bam_caller_vcf_tbi_id_pair,
                                                        truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                        vcfeval_fasta_id, None,
                                                        out_name=bam_caller_out_name if not bed_id else bam_caller_out_name + '-unclipped',
                                                        score_field='GQ', cores=context.config.vcfeval_cores,
                                                        memory=context.config.vcfeval_mem,
                                                        disk=context.config.vcfeval_disk).rv()
                    if happy:
                        happy_results[bam_caller_out_name]["unclipped"] = \
                        bam_caller_job.addFollowOnJobFn(run_happy, context, sample_name, bam_caller_vcf_tbi_id_pair,
                                                        truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                        vcfeval_fasta_id, None,
                                                        out_name=bam_caller_out_name if not bed_id else bam_caller_out_name + '-unclipped',
                                                        cores=context.config.vcfeval_cores,
                                                        memory=context.config.vcfeval_mem,
                                                        disk=context.config.vcfeval_disk).rv()

                    if sveval:
                        sveval_results[bam_caller_out_name]["unclipped"] = \
                        bam_caller_job.addFollowOnJobFn(run_sv_eval, context, sample_name, bam_caller_vcf_tbi_id_pair,
                                                        truth_vcf_id, truth_vcf_tbi_id, min_sv_len=25,
                                                        sv_overlap=0.5, sv_region_overlap=1, bed_id=None,
                                                        out_name=bam_caller_out_name if not bed_id else bam_caller_out_name + '-unclipped',                                                        
                                                        cores=context.config.vcfeval_cores,
                                                        memory=context.config.vcfeval_mem,
                                                        disk=context.config.vcfeval_disk).rv()                        

                vcf_tbi_id_pairs.append(bam_caller_vcf_tbi_id_pair)
                timing_results.append(timing_result)
                names.append(bam_caller_out_name)

    # optional override to filter-opts when running genotype
    # this is allows us to run different filter-opts for call and genotype
    # (something we don't really need an interface for besides in calleval)
    if filter_opts_gt:
        gt_context = copy.deepcopy(context)
        gt_context.config.filter_opts = filter_opts_gt.split()
    else:
        gt_context = context

    if gam_ids:
        for gam_id, gam_idx_id, gam_name, xg_id in zip(gam_ids, gam_idx_ids, gam_names, xg_ids):
            for gt in [False, True]:
                if (call and not gt) or (genotype and gt):
                    out_name = '{}{}'.format(gam_name, '-gt' if gt else '-call')
                    call_job = child_job.addChildJobFn(run_all_calling, gt_context if gt else context,
                                                       xg_id, [gam_id], [gam_idx_id], chroms, vcf_offsets,
                                                       sample_name, genotype=gt, recall=recall,
                                                       out_name=out_name,
                                                       cores=context.config.misc_cores,
                                                       memory=context.config.misc_mem,
                                                       disk=context.config.misc_disk)
                    vcf_tbi_id_pair = (call_job.rv(0), call_job.rv(1))
                    timing_result = call_job.rv(2)

                    if not vcfeval_score_field:
                        score_field = 'GQ' if gt else 'QUAL'
                    else:
                        score_field = vcfeval_score_field

                    if bed_id:
                        eval_results[out_name]["clipped"] = \
                            call_job.addFollowOnJobFn(run_vcfeval, context, sample_name, vcf_tbi_id_pair,
                                                      truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                      vcfeval_fasta_id, bed_id, out_name=out_name,
                                                      score_field=score_field).rv()
                        if happy:
                            happy_results[out_name]["clipped"] = \
                            call_job.addFollowOnJobFn(run_happy, context, sample_name, vcf_tbi_id_pair,
                                                      truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                      vcfeval_fasta_id, bed_id, out_name=out_name).rv()

                        if sveval:
                            sveval_results[out_name]["clipped"] = \
                            call_job.addFollowOnJobFn(run_sv_eval, context, sample_name, vcf_tbi_id_pair,
                                                      truth_vcf_id, truth_vcf_tbi_id, min_sv_len=25,
                                                      sv_overlap=0.5, sv_region_overlap=1,
                                                      bed_id = bed_id, out_name=out_name).rv()
                                                    
                    if not clip_only:
                        # Also do unclipped
                        eval_results[out_name]["unclipped"] = \
                            call_job.addFollowOnJobFn(run_vcfeval, context, sample_name, vcf_tbi_id_pair,
                                                      truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                      vcfeval_fasta_id, None,
                                                      out_name=out_name if not bed_id else out_name + '-unclipped',
                                                      score_field=score_field).rv()
                        if happy:
                            happy_results[out_name]["unclipped"] = \
                            call_job.addFollowOnJobFn(run_happy, context, sample_name, vcf_tbi_id_pair,
                                                      truth_vcf_id, truth_vcf_tbi_id, 'ref.fasta',
                                                      vcfeval_fasta_id, None,
                                                      out_name=out_name if not bed_id else out_name + '-unclipped').rv()

                        if sveval:
                            sveval_results[out_name]["unclipped"] = \
                            call_job.addFollowOnJobFn(run_sv_eval, context, sample_name, vcf_tbi_id_pair,
                                                      truth_vcf_id, truth_vcf_tbi_id, min_sv_len=25,
                                                      sv_overlap=0.5, sv_region_overlap=1,
                                                      bed_id = None,
                                                      out_name=out_name if not bed_id else out_name + '-unclipped').rv()
                            
                    
                    vcf_tbi_id_pairs.append(vcf_tbi_id_pair)
                    timing_results.append(timing_result)
                    names.append(out_name)            
                    

    calleval_results = child_job.addFollowOnJobFn(run_calleval_results, context, names,
                                                  vcf_tbi_id_pairs, eval_results, happy_results, sveval_results,
                                                  timing_results, plot_sets,
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
            inputGamIdxIDs = []
            gamToID = {}
            gaiToID = {}
            if options.gams:
                for gam in options.gams:
                    if gam not in gamToID:
                        gamToID[gam] = toil.importFile(gam)
                    inputGamFileIDs.append(gamToID[gam])
                    gai = gam + '.gai'
                    if gai not in gaiToID:
                        try:
                            gaiToID[gai] = toil.importFile(gai)
                        except:
                            gaiToID[gai] = None
                    inputGamIdxIDs.append(gaiToID[gai])
                        
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
            vcfeval_fasta_id = toil.importFile(options.vcfeval_fasta)
            bed_id = toil.importFile(options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None
            clip_only = options.clip_only
            
            # What do we plot together?
            plot_sets = parse_plot_sets(options.plot_sets) 

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)

            if options.vcfeval_fasta.endswith('.gz'):
                # unzip the fasta for evaluation
                vcfeval_fasta_id = init_job.addChildJobFn(run_unzip_fasta, context, vcfeval_fasta_id,
                                                       os.path.basename(options.vcfeval_fasta)).rv()

            if options.caller_fasta is not None:
                # Calling uses a different FASTA (for a subregion)
                caller_fasta_id = toil.importFile(options.caller_fasta)
                if options.caller_fasta.endswith('.gz'):
                    # unzip the fasta for freebayes
                    caller_fasta_id = init_job.addChildJobFn(run_unzip_fasta, context, caller_fasta_id,
                                                                os.path.basename(options.vcfeval_fasta)).rv()
            else:
                # Use the same FASTA as evaluation
                caller_fasta_id = vcfeval_fasta_id

            if options.chroms:
                # Make sure the comparison respects --chroms if it's provided
                vcf_subset_job = init_job.addChildJobFn(run_vcf_subset, context, vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                                        options.chroms, disk=context.config.vcfeval_disk)
                vcfeval_baseline_id, vcfeval_baseline_tbi_id = vcf_subset_job.rv(0), vcf_subset_job.rv(1)

            # Make a root job
            root_job = Job.wrapJobFn(run_calleval, context, inputXGFileIDs, inputGamFileIDs, inputGamIdxIDs,
                                     inputBamFileIDs, inputBamIdxIds,
                                     options.gam_names, options.bam_names, 
                                     vcfeval_baseline_id, vcfeval_baseline_tbi_id, caller_fasta_id, vcfeval_fasta_id,
                                     bed_id, clip_only,
                                     options.call,
                                     options.genotype,
                                     options.sample_name,
                                     options.chroms, options.vcf_offsets,
                                     options.vcfeval_score_field,
                                     plot_sets,
                                     options.filter_opts_gt,
                                     options.surject,
                                     options.interleaved,
                                     options.freebayes,
                                     options.platypus,
                                     options.happy,
                                     options.sveval,
                                     options.recall,
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
    
    
    
