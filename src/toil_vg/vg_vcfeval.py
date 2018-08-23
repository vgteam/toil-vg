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
    parser.add_argument("--vcfeval", action="store_true",
                        help="run rtg vcfeval comparison.  (will be run by default if no other tool specified)")
    parser.add_argument("--happy", action="store_true",
                        help="run hap.py comparison.")
    parser.add_argument("--sveval", action="store_true",
                        help="run bed-based sv comparison.")
    parser.add_argument("--min_sv_len", type=int, default=20,
                        help="minimum length to consider when doing bed sv comparison (using --sveval)")
    parser.add_argument("--sv_region_overlap", type=float, default=1.0,
                        help="sv must overlap bed region (--vcfeval_bed_regions) by this fraction to be considered")
    parser.add_argument("--sv_overlap", type=float, default=0.5,
                        help="minimum reciprical overlap required for bed intersection to count as TP")
    parser.add_argument("--sv_smooth", type=int, default=0,
                        help="mege up svs (in calls and truth) that are at most this many bases apart")

def validate_vcfeval_options(options):
    """ check some options """
    # we can relax this down the road by optionally doing some compression/indexing
    assert options.vcfeval_baseline and options.vcfeval_baseline.endswith(".vcf.gz")
    assert options.call_vcf.endswith(".vcf.gz")

    assert not (options.happy or options.vcfeval) or options.vcfeval_fasta

    
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

def run_vcfeval_roc_plot(job, context, roc_table_ids, names=[], kind=None, number=0, title=None,
                         show_scores=False, line_width=2, ps_plot=False):
    """
    Draw some rocs from the vcfeval output. Return (snps_id, nonsnps_id,
    weighted_id)
    
    kind specifies the subtype of roc plot (e.g. 'snp-unclipped'). number gives
    the number of this plot in all plots of that kind. title gives an optional
    human-readable title.
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

    # Make sure the kind has 'roc' in it
    if kind is None:
        kind = 'roc'
    else:
        kind = 'roc-{}'.format(kind)
    
    plot_filename = title_to_filename(kind, number, title, 'svg')
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
    """
    
    Run RTG vcf_eval to compare VCFs.
    
    Return a results dict like:
    
    {
        "f1": f1 score as float,
        "summary": summary file ID,
        "archive": output archive ID,
        "snp": ROC .tsv.gz data file ID for SNPs,
        "non_snp": ROC .tsv.gz data file ID for non-SNP variants,
        "weighted": ROC .tsv.gz data file ID for a weighted combination of SNP and non-SNP variants
    }
    
    Some ROC data file IDs may not be present if they were not calculated.
    
    """

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

    # Start the output dict
    out_dict = {
        "f1": f1, 
        "summary": out_summary_id, 
        "archive": out_archive_id
    }

    #  roc data (written to outstore to allow re-plotting)
    for roc_name in ['snp', 'non_snp', 'weighted']:
        roc_file = os.path.join(work_dir, out_tag, '{}_roc.tsv.gz'.format(roc_name))
        if os.path.isfile(roc_file):
            # Save this one
            dest_file = os.path.join('roc', out_tag, '{}_roc.tsv.gz'.format(roc_name))
            out_dict[roc_name] = context.write_output_file(job, roc_file, dest_file)

    return out_dict

def run_happy(job, context, sample, vcf_tbi_id_pair, vcfeval_baseline_id, vcfeval_baseline_tbi_id, 
              fasta_path, fasta_id, bed_id, fasta_idx_id = None, out_name = None):
    """
    
    Run hap.py to compare VCFs.
    
    Return a results dict like:
    
    {
        "parsed_summary": parsed summary file dict, including F1 scores,
        "summary": summary file ID,
        "archive": output archive ID
    }
    
    The parsed summary dict is by variant type ('SNP', 'INDEL'), and then by metric name (like 'METRIC.F1_Score').
    
    """
    

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
    
    return {
        "parsed_summary": happy_results,
        "summary": out_summary_id,
        "archive": out_archive_id
    }

def run_sv_eval(job, context, sample, vcf_tbi_id_pair, vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                min_sv_len, sv_overlap, sv_region_overlap, sv_smooth = 0, bed_id = None,  out_name = ''):
    """ Run something like Peter Audano's bed-based comparison.  Uses bedtools and bedops to do
    overlap comparison between indels.  Of note: the actual sequence of insertions is never checked!"""

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download the vcf
    call_vcf_id, call_tbi_id = vcf_tbi_id_pair[0], vcf_tbi_id_pair[1]
    call_vcf_name = "{}calls.vcf.gz".format(out_name)
    job.fileStore.readGlobalFile(vcf_tbi_id_pair[0], os.path.join(work_dir, call_vcf_name))
    job.fileStore.readGlobalFile(vcf_tbi_id_pair[1], os.path.join(work_dir, call_vcf_name + '.tbi'))

    # and the truth vcf
    vcfeval_baseline_name = '{}truth.vcf.gz'.format(out_name)
    job.fileStore.readGlobalFile(vcfeval_baseline_id, os.path.join(work_dir, vcfeval_baseline_name))
    job.fileStore.readGlobalFile(vcfeval_baseline_tbi_id, os.path.join(work_dir, vcfeval_baseline_name + '.tbi'))

    # and the regions bed
    if bed_id:
        regions_bed_name = '{}regions.bed'.format(out_name)
        job.fileStore.readGlobalFile(bed_id, os.path.join(work_dir, regions_bed_name))
    else:
        regions_bed_name = None

    if out_name and not out_name.endswith('_'):
        out_name = '{}_'.format(out_name)

    # convert vcfs to BEDs
    call_bed_name = '{}calls.bed'.format(out_name)
    baseline_bed_name = '{}truth.bed'.format(out_name)
    for vcf_name, bed_name in zip([call_vcf_name, vcfeval_baseline_name], [call_bed_name, baseline_bed_name]):
        # vg call sometimes makes IDs that are too long for bedops.  we zap the ids here just in case
        temp_vcf_name = 'idfix-{}.vcf'.format(os.path.splitext(os.path.basename(vcf_name))[0])
        with open(os.path.join(work_dir, temp_vcf_name), 'w') as temp_vcf_file:
            context.runner.call(job, ['bcftools', 'view', '-h', vcf_name],
                                work_dir = work_dir, outfile = temp_vcf_file)
            context.runner.call(job, [['bcftools', 'view', '-H', vcf_name],
                                      ['awk', '-F', '\t', '-v', 'OFS=\t', '{$3=\".\"; print $0;}']],
                                work_dir = work_dir, outfile = temp_vcf_file)
        # then convert the fixed if vcf into bed with bedops
        with open(os.path.join(work_dir, bed_name + '.1del'), 'w') as bed_file:
            vcf2bed_cmd = [['cat', temp_vcf_name], ['vcf2bed']]
            context.runner.call(job, vcf2bed_cmd, work_dir = work_dir, outfile = bed_file, tool_name = 'bedops')
            
        # then expand the deletions, as vcf2bed writes them all with 1-sized intervals for some reason
        expand_deletions(os.path.join(work_dir, bed_name + '.1del'), os.path.join(work_dir, bed_name))

    # if the region bed's there, intersect both the calls and baseline       
    if bed_id:
        clipped_calls_name = "{}clipped_calls.bed".format(out_name)
        clipped_baseline_name = "{}clipped_baseline.bed".format(out_name)

        """
        bedtools intersect -a SV_FILE.bed -b sd_regions_200_0.bed -wa -f 0.50 -u

        * Get all SV variants that intersect SDs
        * -wa: Print whole SV record (not just the intersecting part)
        * -f 0.50: Require 50% of the SV to be in the SD
        * Prevents very large DELs from being tagged if they intersect with part of a region.
        * -u: Print matching SV call only once (even if it intersects with multiple variants)
        * Probably has no effect because of -f 0.50 and the way merging is done, but I leave it in in case I tweak something.
        """
        for bed_name, clipped_bed_name in zip([call_bed_name, baseline_bed_name],
                                              [clipped_calls_name, clipped_baseline_name]):
            with open(os.path.join(work_dir, clipped_bed_name), 'w') as clipped_bed_file:
                context.runner.call(job, ['bedtools', 'intersect', '-a', bed_name, '-b', regions_bed_name,
                                          '-wa', '-f', str(sv_region_overlap), '-u'],
                                    work_dir = work_dir, outfile = clipped_bed_file)
    else:
        clipped_calls_name = call_bed_name
        clipped_baseline_name = baseline_bed_name
        
    # now do the intersection comparison

    """
    To compare SV calls among sets, I typically use a 50% reciprocal overlap. A BED record of an insertion is only the point of insertion (1 bp) regardless of how large the insertion is. To compare insertions, I add the SV length (SVLEN) to the start position and use bedtools overlap. WARNING: Remember to collapse the insertion records back to one bp. If you accidentally intersected with SDs or TRF regions, SVs would arbitrarily hit records they don't actually intersect. I have snakemake pipelines handle this (bed files in "byref" or "bylen" directories) so it can never be mixed up.

    Intersect SVs:

    bedtools intersect -a SV_FILE_A.bed -b SV_FILE_B -wa -f 0.50 -r -u

    * Same bedtools command as before, but with "-r" to force A to overlap with B by 50% AND B to overlap with A by 50%. "-u" becomes important because clustered insertions (common in tandem repeats) will overlap with more than one variant.
    """

    # expand the bed regions to include the insertion lengths for both sets.
    # we also break insertions and deletions into separate file.  Otherwise, insertions can
    # intersect with deletions, inflating our accuracy.
    clipped_calls_ins_name = '{}ins_calls.bed'.format(out_name)
    clipped_baseline_ins_name = '{}ins_calls_baseline.bed'.format(out_name)
    clipped_calls_del_name = '{}del_calls.bed'.format(out_name)
    clipped_baseline_del_name = '{}del_calls_baseline.bed'.format(out_name)
    expand_insertions(os.path.join(work_dir, clipped_calls_name), os.path.join(work_dir, clipped_calls_ins_name),
                      os.path.join(work_dir, clipped_calls_del_name), min_sv_len)
    expand_insertions(os.path.join(work_dir, clipped_baseline_name), os.path.join(work_dir, clipped_baseline_ins_name),
                      os.path.join(work_dir, clipped_baseline_del_name), min_sv_len)

    tp_ins_name = '{}ins-TP-call.bed'.format(out_name)
    tp_del_name = '{}del-TP-call.bed'.format(out_name)
    tp_ins_rev_name = '{}ins-TP-baseline.bed'.format(out_name)
    tp_del_rev_name = '{}del-TP-baseline.bed'.format(out_name)
    fp_ins_name = '{}ins-FP.bed'.format(out_name)
    fp_del_name = '{}del-FP.bed'.format(out_name)
    fn_ins_name = '{}ins-FN.bed'.format(out_name)
    fn_del_name = '{}del-FN.bed'.format(out_name)        
    for tp_indel_name, tp_indel_rev_name, calls_bed_indel_name, baseline_bed_indel_name, fp_indel_name, fn_indel_name in \
        [(tp_ins_name, tp_ins_rev_name, clipped_calls_ins_name, clipped_baseline_ins_name, fp_ins_name, fn_ins_name),
         (tp_del_name, tp_del_rev_name, clipped_calls_del_name, clipped_baseline_del_name, fp_del_name, fn_del_name)]:

        # smooth out features so, say, two side-by-side deltions get treated as one.  this is in keeping with
        # the coarse-grained nature of the analysis and is optional
        if sv_smooth > 0:
            calls_merge_name = calls_bed_indel_name[:-4] + '_merge.bed'
            baseline_merge_name = baseline_bed_indel_name[:-4] + '_merge.bed'            
            with open(os.path.join(work_dir, calls_merge_name), 'w') as calls_merge_file:
                context.runner.call(job, ['bedtools', 'merge', '-d', str(sv_smooth), '-i', calls_bed_indel_name],
                                    work_dir = work_dir, outfile = calls_merge_file)
            with open(os.path.join(work_dir, baseline_merge_name), 'w') as baseline_merge_file:
                context.runner.call(job, ['bedtools', 'merge', '-d', str(sv_smooth), '-i', baseline_bed_indel_name],
                                    work_dir = work_dir, outfile = baseline_merge_file)
            calls_bed_indel_name = calls_merge_name
            baseline_bed_indel_name = baseline_merge_name
        
        # run the 50% overlap test described above for insertions and deletions
        with open(os.path.join(work_dir, tp_indel_name), 'w') as tp_file:
            context.runner.call(job, ['bedtools', 'intersect', '-a', calls_bed_indel_name,
                                  '-b', baseline_bed_indel_name, '-wa', '-f', str(sv_overlap), '-r', '-u'],
                            work_dir = work_dir, outfile = tp_file)            
        # we run other way so we can get false positives from the calls and true positives from the baseline    
        with open(os.path.join(work_dir, tp_indel_rev_name), 'w') as tp_file:        
            context.runner.call(job, ['bedtools', 'intersect', '-b', calls_bed_indel_name,
                                  '-a', baseline_bed_indel_name, '-wa', '-f', str(sv_overlap), '-r', '-u'],
                            work_dir = work_dir, outfile = tp_file)
        # put the false positives in their own file
        with open(os.path.join(work_dir, fp_indel_name), 'w') as fp_file:
            context.runner.call(job, ['bedtools', 'subtract', '-a', calls_bed_indel_name,
                                      '-b', tp_indel_name], work_dir = work_dir, outfile = fp_file)
        # and the false negatives
        with open(os.path.join(work_dir, fn_indel_name), 'w') as fn_file:
            context.runner.call(job, ['bedtools', 'subtract', '-a', baseline_bed_indel_name,
                                      '-b', tp_indel_rev_name], work_dir = work_dir, outfile = fn_file)
        # todo: should we write them out in vcf as well?

    # summarize results into a table
    results = summarize_sv_results(os.path.join(work_dir, tp_ins_name),
                                   os.path.join(work_dir, tp_ins_rev_name),
                                   os.path.join(work_dir, fp_ins_name),
                                   os.path.join(work_dir, fn_ins_name),
                                   os.path.join(work_dir, tp_del_name),
                                   os.path.join(work_dir, tp_del_rev_name),
                                   os.path.join(work_dir, fp_del_name),
                                   os.path.join(work_dir, fn_del_name))

    # write the results to a file
    summary_name = os.path.join(work_dir, '{}sv_accuracy.tsv'.format(out_name))
    with open(summary_name, 'w') as summary_file:
        header = ['Cat', 'TP', 'TP-baseline', 'FP', 'FN', 'Precision', 'Recall', 'F1']
        summary_file.write('#' + '\t'.join(header) +'\n')
        summary_file.write('\t'.join(str(x) for x in ['Total'] + [results[c] for c in header[1:]]) + '\n')
        summary_file.write('\t'.join(str(x) for x in ['INS'] + [results['{}-INS'.format(c)] for c in header[1:]]) + '\n')
        summary_file.write('\t'.join(str(x) for x in ['DEL'] + [results['{}-DEL'.format(c)] for c in header[1:]]) + '\n')
    summary_id = context.write_output_file(job, os.path.join(work_dir, summary_name))

    # tar up some relavant data
    tar_dir = os.path.join(work_dir, '{}sv_evaluation'.format(out_name))
    os.makedirs(tar_dir)
    for dir_file in os.listdir(work_dir):
        if os.path.splitext(dir_file)[1] in ['.bed', '.tsv', '.vcf.gz', '.vcf.gz.tbi']:
            shutil.copy2(os.path.join(work_dir, dir_file), os.path.join(tar_dir, dir_file))
    context.runner.call(job, ['tar', 'czf', os.path.basename(tar_dir) + '.tar.gz', os.path.basename(tar_dir)],
                        work_dir = work_dir)
    archive_id = context.write_output_file(job, os.path.join(work_dir, tar_dir + '.tar.gz'))

    return results

def expand_deletions(in_bed_name, out_bed_name):
    """
    Expand every deletion in a BED file and update its end coordinate to reflect its length.
    with open(in_bed_name) as in_bed, open(out_ins_bed_name, 'w') as out_ins, open(out_del_bed_name, 'w') as out_del:
    ** We do this before intersection with regions of interest **
    """
    with open(in_bed_name) as in_bed, open(out_bed_name, 'w') as out_bed:
        for line in in_bed:
            if line.strip():
                toks = line.strip().split('\t')
                ref_len = len(toks[5])
                alt_len = len(toks[6])
                if ref_len > alt_len:
                    assert int(toks[2]) == int(toks[1]) + 1
                    # expand the deletion
                    toks[2] = str(int(toks[1]) + ref_len)
                    out_bed.write('\t'.join(toks) + '\n')
                else:
                    # leave insertions as is for now
                    out_bed.write(line)

def expand_insertions(in_bed_name, out_ins_bed_name, out_del_bed_name, min_sv_size):
    """
    Go through every insertion in a BED file and update its end coordinate to reflect its length.  This is done
    to compare two sets of insertions.  We also break out deletions into their own file. 
    ** We do this after intersection with regions of interest but before comparison **
    """
    with open(in_bed_name) as in_bed, open(out_ins_bed_name, 'w') as out_ins, open(out_del_bed_name, 'w') as out_del:
        for line in in_bed:
            if line.strip():
                toks = line.strip().split('\t')
                ref_len = len(toks[5])
                alt_len = len(toks[6])
                # filter out some vg call nonsense
                if toks[5] == '.' or toks[6] == '.':
                    continue
                if ref_len < alt_len and alt_len >= min_sv_size:
                    assert int(toks[2]) == int(toks[1]) + 1
                    # expand the insertion
                    toks[2] = str(int(toks[1]) + alt_len)
                    out_ins.write('\t'.join(toks) + '\n')
                elif ref_len >= min_sv_size:
                    # just filter out the deletion
                    out_del.write(line)

def summarize_sv_results(tp_ins, tp_ins_baseline, fp_ins, fn_ins,
                         tp_del, tp_del_baseline, fp_del, fn_del):
    """
    Use the various bed files to compute accuracies.  Also return a tarball of 
    all the files used.
    """
    def wc(f):
        with open(f) as ff:
            return sum(1 for line in ff)

    def pr(tp, fp, fn):
        prec = float(tp) / float(tp + fp) if tp + fp else 0
        rec = float(tp) / float(tp + fn) if tp + fn else 0
        f1 = 2.0 * tp / float(2 * tp + fp + fn) if tp else 0
        return prec, rec, f1

    header = []
    row = []

    # results in dict
    results = {}
    results['TP-INS'] = wc(tp_ins)
    results['TP-baseline-INS'] = wc(tp_ins_baseline)
    results['FP-INS'] = wc(fp_ins)
    results['FN-INS'] = wc(fn_ins)
    ins_pr = pr(results['TP-INS'], results['FP-INS'], results['FN-INS'])
    results['Precision-INS'] = ins_pr[0]
    results['Recall-INS'] = ins_pr[1]
    results['F1-INS'] = ins_pr[2]

    results['TP-DEL'] = wc(tp_del)
    results['TP-baseline-DEL'] = wc(tp_del_baseline)
    results['FP-DEL'] = wc(fp_del)
    results['FN-DEL'] = wc(fn_del)
    del_pr = pr(results['TP-DEL'], results['FP-DEL'], results['FN-DEL'])
    results['Precision-DEL'] = del_pr[0]
    results['Recall-DEL'] = del_pr[1]
    results['F1-DEL'] = del_pr[2]

    results['TP'] = results['TP-INS'] + results['TP-DEL']
    results['TP-baseline'] = results['TP-baseline-INS'] + results['TP-baseline-DEL']    
    results['FP'] = results['FP-INS'] + results['FP-DEL']
    results['FN'] = results['FN-INS'] + results['FN-DEL']
    tot_pr = pr(results['TP'], results['FP'], results['FN'])
    results['Precision'] = tot_pr[0]
    results['Recall'] = tot_pr[1]
    results['F1'] = tot_pr[2]

    return results                    

def vcfeval_main(context, options):
    """ command line access to toil vcf eval logic"""

    # check some options
    validate_vcfeval_options(options)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    # Default to vcfeval
    if not options.happy and not options.sveval:
        options.vcfeval = True

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:
            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            vcfeval_baseline_id = toil.importFile(options.vcfeval_baseline)
            call_vcf_id = toil.importFile(options.call_vcf)
            vcfeval_baseline_tbi_id = toil.importFile(options.vcfeval_baseline + '.tbi')
            call_tbi_id = toil.importFile(options.call_vcf + '.tbi')            
            fasta_id = toil.importFile(options.vcfeval_fasta) if options.vcfeval_fasta else None
            bed_id = toil.importFile(options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)
            
            # Make a root job
            if options.vcfeval:
                vcfeval_job = Job.wrapJobFn(run_vcfeval, context, None,
                                            (call_vcf_id, call_tbi_id),
                                            vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                            options.vcfeval_fasta, fasta_id, bed_id,
                                            score_field=options.vcfeval_score_field,
                                            cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                            disk=context.config.vcfeval_disk)
                init_job.addFollowOn(vcfeval_job)
                
            if options.happy:
                happy_job = Job.wrapJobFn(run_happy, context, None,
                                          (call_vcf_id, call_tbi_id),
                                          vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                          options.vcfeval_fasta, fasta_id, bed_id,
                                          cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                          disk=context.config.vcfeval_disk)
                init_job.addFollowOn(happy_job)

            if options.sveval:                
                sv_job = Job.wrapJobFn(run_sv_eval, context, None,
                                       (call_vcf_id, call_tbi_id),
                                       vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                       options.min_sv_len, options.sv_overlap, options.sv_region_overlap,
                                       options.sv_smooth, bed_id, 
                                       cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                       disk=context.config.vcfeval_disk)
                init_job.addFollowOn(sv_job)

            # Run the job
            toil.start(init_job)
        else:
            toil.restart()

    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
