#!/usr/bin/env python2.7
"""
Thin wrapper for vcfeval, as a convenience to stick a vcfeval output directory
along with the other toil-vg output.  Can be run standalone as well.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, time, timeit, errno
from uuid import uuid4
import logging

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore

logger = logging.getLogger(__name__)

def vcfeval_subparser(parser):
    """
    Create a subparser for vcf evaluation.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("--call_vcf", type=make_url, required=True,
                        help="input vcf (must be bgzipped and have .tbi")
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common calling options shared with toil_vg pipeline
    vcfeval_parse_args(parser)

    # Add common docker options shared with toil_vg pipeline
    add_container_tool_parse_args(parser)

def vcfeval_parse_args(parser):
    """ centralize reusable vcfevaling parameters here """

    parser.add_argument("--vcfeval_baseline", type=make_url,
                        help="Path to baseline VCF file for comparison (must be bgzipped and have .tbi)")

    parser.add_argument("--vcfeval_fasta", type=make_url,
                        help="Path to DNA sequence file, required for vcfeval. Maybe be gzipped")

    parser.add_argument("--vcfeval_bed_regions", type=make_url,
                        help="BED file of regions to consider")
    parser.add_argument("--vcfeval_opts", type=str,
                        help="Additional options for vcfeval (wrapped in \"\")",
                        default=None)
    parser.add_argument("--vcfeval_cores", type=int,
                        default=1,
                        help="Cores to use for vcfeval")
    parser.add_argument("--vcfeval_score_field", default=None,
                        help="vcf FORMAT field to use for ROC score.  overrides vcfeval_opts")
    parser.add_argument("--happy", action="store_true",
                        help="run hap.py comparison in addition to rtg vcfeval")

def validate_vcfeval_options(options):
    """ check some options """
    # we can relax this down the road by optionally doing some compression/indexing
    assert options.vcfeval_baseline and options.vcfeval_baseline.endswith(".vcf.gz")
    assert options.call_vcf.endswith(".vcf.gz")

    assert options.vcfeval_fasta

    
def parse_f1(summary_path):
    """ grab the best f1 out of vcfeval's summary.txt """

    # get the F1 out of summary.txt
    # expect header on 1st line and data on 3rd and below
    # we take the best F1 found over these lines (which should correspond to best
    # point on quality ROC curve)
    # todo: be more robust
    f1 = None
    with open(summary_path) as sum_file:
        header = sum_file.readline().split()
        assert header[-1] == 'F-measure'        
        line = sum_file.readline()
        for line in sum_file:
            data = line.strip().split()
            assert len(data) == len(header)
            line_f1 = float(data[-1])
            if f1 is None or line_f1 > f1:
                f1 = line_f1
    return f1

def parse_happy_summary(summary_path):
    """ Turn the happy summary into a dictionary we can easily get F1 scores etc. from """
    results = {}  # make something like: results[INDEL][METRIC.F1_Score] = 0.9x
    with open(summary_path) as sum_file:
        header = sum_file.readline().split(',')
        for line in sum_file:
            row = line.split(',')
            cat = row[0]
            if row[1] == 'ALL':
                cat += '.ALL'
            assert cat not in results
            results[cat] = {}
            for column in range(1, len(header)):
                results[cat][header[column]] = row[column] if len(row[column]) else '0'
        return results

def run_vcfeval_roc_plot(job, context, roc_table_ids, names=[], title=None, show_scores=False,
                          line_width=2, ps_plot=False):
    """
    draw some rocs from the vcfeval output. return (snps_id, nonsnps_id, weighted_id)
    """

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # dummy default names
    if not names:
        names = ['vcfeval_output-{}'.format(i) for i in range(len(roc_table_ids))]

    # rely on unique input names
    assert len(names) == len(set(names))
        
    # download the files
    table_file_paths = [os.path.join(work_dir, name, name) + '.tsv.gz' for name in names]
    table_file_rel_paths = [os.path.join(name, name) + '.tsv.gz' for name in names]
    for table_path, name, file_id in zip(table_file_paths, names, roc_table_ids):
        # rtg gets naming information from directory structure, so we read each file into
        # its own dir
        os.makedirs(os.path.join(work_dir, name))
        job.fileStore.readGlobalFile(file_id, table_path)
        RealtimeLogger.info('Downloaded {} to {}'.format(file_id, table_path))

    plot_filename = 'roc{}.svg'.format('-{}'.format(title) if title else '')
    out_roc_path = os.path.join(work_dir, plot_filename)

    roc_opts = []
    if title:
        roc_opts += ['--title', title]
    if show_scores:
        roc_opts += ['--scores']
    if line_width:
        roc_opts += ['--line-width', line_width]
    if ps_plot:
        roc_opts += ['precision-sensitivity']

    out_ids = []

    roc_cmd = ['rtg', 'rocplot', '--svg', os.path.basename(out_roc_path)]
    roc_cmd += roc_opts + table_file_rel_paths
    
    context.runner.call(job, roc_cmd, work_dir = work_dir)

    return context.write_output_file(job, out_roc_path, os.path.join('plots', plot_filename))

def run_extract_sample_truth_vcf(job, context, sample, input_baseline_id, input_baseline_tbi_id):
    """
    
    Extract a single-sample truth VCF from the given truth VCF .vcf.gz and .vcf.gz.tbi.
    
    Returns a pair of file IDs for the resulting .vcf.gz and .vcf.gz.tbi.
    
    Filtering the truth down is useful because it can save memory.
    
    TODO: use this in toil-vg vcfeval, instead of just providing it as a utility.
    
    """
    
    # Make a local work directory
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the truth vcf
    vcfeval_baseline_name = 'full-truth.vcf.gz'
    job.fileStore.readGlobalFile(input_baseline_id, os.path.join(work_dir, vcfeval_baseline_name))
    job.fileStore.readGlobalFile(input_baseline_tbi_id, os.path.join(work_dir, vcfeval_baseline_name + '.tbi'))    
    
    # Make the single-sample VCF
    single_name = 'single-truth-{}.vcf.gz'.format(sample)
    context.runner.call(job, ['bcftools', 'view', vcfeval_baseline_name, '-s', sample, '-o', single_name, '-O', 'z'], work_dir=work_dir)
    
    # Index it
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', single_name], work_dir=work_dir)
    
    # Upload file and index
    single_vcf_id = context.write_output_file(job, os.path.join(work_dir, single_name))
    single_vcf_tbi_id = context.write_output_file(job, os.path.join(work_dir, single_name + '.tbi'))
    
    return single_vcf_id, single_vcf_tbi_id

def run_vcfeval(job, context, sample, vcf_tbi_id_pair, vcfeval_baseline_id, vcfeval_baseline_tbi_id, 
    fasta_path, fasta_id, bed_id, out_name = None, score_field=None):
    """ run vcf_eval, return (f1 score, summary id, output archive id, snp-id, nonsnp-id, weighted-id)"""

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download the vcf
    call_vcf_id, call_tbi_id = vcf_tbi_id_pair[0], vcf_tbi_id_pair[1]
    call_vcf_name = "calls.vcf.gz"
    job.fileStore.readGlobalFile(vcf_tbi_id_pair[0], os.path.join(work_dir, call_vcf_name))
    job.fileStore.readGlobalFile(vcf_tbi_id_pair[1], os.path.join(work_dir, call_vcf_name + '.tbi'))

    # and the truth vcf
    vcfeval_baseline_name = 'truth.vcf.gz'
    job.fileStore.readGlobalFile(vcfeval_baseline_id, os.path.join(work_dir, vcfeval_baseline_name))
    job.fileStore.readGlobalFile(vcfeval_baseline_tbi_id, os.path.join(work_dir, vcfeval_baseline_name + '.tbi'))    
    # download the fasta (make sure to keep input extension)
    fasta_name = "fa_" + os.path.basename(fasta_path)
    job.fileStore.readGlobalFile(fasta_id, os.path.join(work_dir, fasta_name))

    # download the bed regions
    bed_name = "bed_regions.bed" if bed_id else None
    if bed_id:
        job.fileStore.readGlobalFile(bed_id, os.path.join(work_dir, bed_name))

    # use out_name if specified, otherwise sample
    if sample and not out_name:
        out_name = sample        
    if out_name:
        out_tag = '{}_vcfeval_output'.format(out_name)
    else:
        out_tag = 'vcfeval_output'
        
    # output directory
    out_name = out_tag
    # indexed sequence
    sdf_name = fasta_name + ".sdf"
    
    # make an indexed sequence (todo: allow user to pass one in)
    context.runner.call(job, ['rtg', 'format',  fasta_name, '-o', sdf_name], work_dir=work_dir)    

    # run the vcf_eval command
    cmd = ['rtg', 'vcfeval', '--calls', call_vcf_name,
           '--baseline', vcfeval_baseline_name,
           '--template', sdf_name, '--output', out_name,
           '--threads', str(context.config.vcfeval_cores)]

    if bed_name is not None:
        cmd += ['--evaluation-regions', bed_name]

    if context.config.vcfeval_opts:
        cmd += context.config.vcfeval_opts

    # override score field from options with one from parameter
    if score_field:
        for opt in ['-f', '--vcf-score-field']:
            if opt in cmd:
                opt_idx = cmd.index(opt)
                del cmd[opt_idx]
                del cmd[opt_idx]
        cmd += ['--vcf-score-field', score_field]
        
    if sample:
        # Pass the sample name along, since it is needed if the truth VCF has multiple samples
        cmd += ['--sample', sample]

    
    try:
        context.runner.call(job, cmd, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the alignment
        logging.error("VCF evaluation failed. Dumping files.")
        context.write_output_file(job, os.path.join(work_dir, call_vcf_name))
        context.write_output_file(job, os.path.join(work_dir, vcfeval_baseline_name))
        # TODO: Dumping the sdf folder doesn't seem to work right. But we can dump the fasta
        context.write_output_file(job, os.path.join(work_dir, fasta_name))
        if bed_name is not None:
            context.write_output_file(job, os.path.join(work_dir, bed_name))
        
        raise


    # copy results to outstore 
    
    # vcfeval_output_summary.txt
    out_summary_id = context.write_output_file(job, os.path.join(work_dir, out_tag, 'summary.txt'),
                                               out_store_path = '{}_summary.txt'.format(out_tag))

    # vcfeval_output.tar.gz -- whole shebang
    context.runner.call(job, ['tar', 'czf', out_tag + '.tar.gz', out_tag], work_dir = work_dir)
    out_archive_id = context.write_output_file(job, os.path.join(work_dir, out_tag + '.tar.gz'))

    # truth VCF
    context.write_output_file(job, os.path.join(work_dir, vcfeval_baseline_name))
    context.write_output_file(job, os.path.join(work_dir, vcfeval_baseline_name + '.tbi'))
    
    # vcfeval_output_f1.txt (used currently by tests script)
    f1 = parse_f1(os.path.join(work_dir, os.path.basename(out_name), "summary.txt"))
    f1_path = os.path.join(work_dir, "f1.txt")    
    with open(f1_path, "w") as f:
        f.write(str(f1))
    context.write_output_file(job, f1_path, out_store_path = '{}_f1.txt'.format(out_tag))

    #  roc data (written to outstore to allow re-plotting)
    out_roc_ids = []
    for roc_name in ['snp', 'non_snp', 'weighted']:
        roc_file = os.path.join(work_dir, out_tag, '{}_roc.tsv.gz'.format(roc_name))
        if os.path.isfile(roc_file):
            out_roc_ids.append(context.write_output_file(job, roc_file,
                               os.path.join('roc', out_tag, '{}_roc.tsv.gz'.format(roc_name))))
        else:
            out_roc_ids.append(None)

    return [f1, out_summary_id, out_archive_id] + out_roc_ids

def run_happy(job, context, sample, vcf_tbi_id_pair, vcfeval_baseline_id, vcfeval_baseline_tbi_id, 
              fasta_path, fasta_id, bed_id, fasta_idx_id = None, out_name = None):
    """ run vcf_eval, return (f1 score, summary id, output archive id, snp-id, nonsnp-id, weighted-id)"""

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download the vcf
    call_vcf_id, call_tbi_id = vcf_tbi_id_pair[0], vcf_tbi_id_pair[1]
    call_vcf_name = "calls.vcf.gz"
    job.fileStore.readGlobalFile(vcf_tbi_id_pair[0], os.path.join(work_dir, call_vcf_name))
    job.fileStore.readGlobalFile(vcf_tbi_id_pair[1], os.path.join(work_dir, call_vcf_name + '.tbi'))

    # and the truth vcf
    vcfeval_baseline_name = 'truth.vcf.gz'
    job.fileStore.readGlobalFile(vcfeval_baseline_id, os.path.join(work_dir, vcfeval_baseline_name))
    job.fileStore.readGlobalFile(vcfeval_baseline_tbi_id, os.path.join(work_dir, vcfeval_baseline_name + '.tbi'))
    
    # download the fasta (make sure to keep input extension)
    fasta_name = "fa_" + os.path.basename(fasta_path)
    job.fileStore.readGlobalFile(fasta_id, os.path.join(work_dir, fasta_name))

    # download or create the fasta index which is required by hap.py
    if fasta_idx_id:
        job.fileStore.readGlobalFile(fasta_idx_id, os.path.join(work_dir, fasta_name + '.fai'))
    else:
        context.runner.call(job, ['samtools', 'faidx', fasta_name], work_dir=work_dir)

    # download the bed regions
    bed_name = "bed_regions.bed" if bed_id else None
    if bed_id:
        job.fileStore.readGlobalFile(bed_id, os.path.join(work_dir, bed_name))

    # use out_name if specified, otherwise sample
    if sample and not out_name:
        out_name = sample        
    if out_name:
        out_tag = '{}_happy'.format(out_name)
    else:
        out_tag = 'happy_output'
        
    # output directory
    out_name = os.path.join(out_tag, 'happy')
    os.makedirs(os.path.join(work_dir, out_tag))
    
    # run the hap.py command
    cmd = ['hap.py', vcfeval_baseline_name, call_vcf_name,
           '--report-prefix', out_name, '--reference', fasta_name, '--write-vcf', '--write-counts', '--no-roc',
           '--threads', str(job.cores)]

    if bed_name:
        cmd += ['--false-positives', bed_name]
        
    try:
        context.runner.call(job, cmd, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the alignment
        logging.error("hap.py VCF evaluation failed. Dumping files.")
        context.write_output_file(job, os.path.join(work_dir, call_vcf_name))
        context.write_output_file(job, os.path.join(work_dir, vcfeval_baseline_name))
        # TODO: Dumping the sdf folder doesn't seem to work right. But we can dump the fasta
        context.write_output_file(job, os.path.join(work_dir, fasta_name))
        if bed_name is not None:
            context.write_output_file(job, os.path.join(work_dir, bed_name))
        
        raise


    # copy results to outstore 
    
    # happy_output_summary.csv
    out_summary_id = context.write_output_file(job, os.path.join(work_dir, out_name + '.summary.csv'),
                                               out_store_path = '{}_summary.csv'.format(out_tag))

    # happy_output.tar.gz -- whole shebang
    context.runner.call(job, ['tar', 'czf', out_tag + '.tar.gz', out_tag], work_dir = work_dir)
    out_archive_id = context.write_output_file(job, os.path.join(work_dir, out_tag + '.tar.gz'))

    # happy_output_f1.txt Snp1-F1 TAB Indel-F1 (of variants marked as PASS)
    happy_results = parse_happy_summary(os.path.join(work_dir, out_name + '.summary.csv'))
    f1_path = os.path.join(work_dir, "happy_f1.txt")    
    with open(f1_path, "w") as f:
        f.write('{}\t{}\n'.format(happy_results['SNP']['METRIC.F1_Score'], happy_results['INDEL']['METRIC.F1_Score']))
    context.write_output_file(job, f1_path, out_store_path = '{}_f1.txt'.format(out_tag))
    
    return happy_results, out_summary_id, out_archive_id

def vcfeval_main(context, options):
    """ command line access to toil vcf eval logic"""

    # check some options
    validate_vcfeval_options(options)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:
            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            vcfeval_baseline_id = toil.importFile(options.vcfeval_baseline)
            call_vcf_id = toil.importFile(options.call_vcf)
            vcfeval_baseline_tbi_id = toil.importFile(options.vcfeval_baseline + '.tbi')
            call_tbi_id = toil.importFile(options.call_vcf + '.tbi')            
            fasta_id = toil.importFile(options.vcfeval_fasta)
            bed_id = toil.importFile(options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            vcfeval_job = Job.wrapJobFn(run_vcfeval, context, None,
                                        (call_vcf_id, call_tbi_id),
                                        vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                        options.vcfeval_fasta, fasta_id, bed_id,
                                        score_field=options.vcfeval_score_field,
                                        cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                        disk=context.config.vcfeval_disk)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)
            init_job.addFollowOn(vcfeval_job)
            if options.happy:
                happy_job = Job.wrapJobFn(run_happy, context, None,
                                          (call_vcf_id, call_tbi_id),
                                          vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                          options.vcfeval_fasta, fasta_id, bed_id,
                                          cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                          disk=context.config.vcfeval_disk)
                init_job.addFollowOn(happy_job)

            # Run the job
            f1 = toil.start(init_job)
        else:
            f1 = toil.restart()

        print("F1 Score : {}".format(f1))
                
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
