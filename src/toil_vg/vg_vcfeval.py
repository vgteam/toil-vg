#!/usr/bin/env python2.7
"""
Thin wrapper for vcfeval, as a convenience to stick a vcfeval output directory
along with the other toil-vg output.  Can be run standalone as well.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, time, timeit, errno, vcf
from uuid import uuid4
import logging

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_call import sort_vcf
from toil_vg.vg_construct import run_make_control_vcfs

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
    parser.add_argument("--max_sv_len", type=int, default=sys.maxint,
                        help="maximum length to consider when doing bed sv comparison (using --sveval)")    
    parser.add_argument("--sv_region_overlap", type=float, default=1.0,
                        help="sv must overlap bed region (--vcfeval_bed_regions) by this fraction to be considered")
    parser.add_argument("--sv_overlap", type=float, default=0.5,
                        help="minimum overlap coverage required for bed intersection to count as TP")
    parser.add_argument("--ins_max_gap", type=int, default=20,
                        help="maximum distance between insertions to be compared")
    parser.add_argument("--ins_seq_comp", action="store_true",
                        help="compare insertion sequence instead of their size only.")
    parser.add_argument("--del_min_rol", type=float, default=0.1,
                        help="the minimum reciprocal overlap when computing coverage on deletions")
    parser.add_argument("--check_inv", action="store_true",
                        help="try to identify inversions during SV evaluation.")
    parser.add_argument("--genotype_eval", action="store_true",
                        help="genotype evaluation instead of calling evaluation during SV evaluation.")
    parser.add_argument("--normalize", action="store_true",
                        help="normalize both VCFs before SV comparison with bcftools norm (requires --vcfeva_fasta)")
    parser.add_argument("--normalize_baseline", action="store_true",
                        help="normalize the truth VCF before SV comparison with bcftools norm (requires --vcfeva_fasta)")
    parser.add_argument("--normalize_calls", action="store_true",
                        help="normalize the calls VCFs before SV comparison with bcftools norm (requires --vcfeva_fasta)")
    parser.add_argument("--vcfeval_sample",
                        help="extract this sample from calls and truth vcf (if possible) before comparison")
                        

def validate_vcfeval_options(options):
    """ check some options """
    # we can relax this down the road by optionally doing some compression/indexing
    assert options.vcfeval_baseline and options.vcfeval_baseline.endswith(".vcf.gz")
    assert options.call_vcf.endswith(".vcf.gz")

    assert not (options.happy or options.vcfeval) or options.vcfeval_fasta
    assert not options.normalize or options.vcfeval_fasta
    assert not options.normalize_calls or options.vcfeval_fasta
    assert not options.normalize_baseline or options.vcfeval_fasta

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
                min_sv_len, max_sv_len, sv_overlap, sv_region_overlap, bed_id = None,
                ins_ref_len=10, del_min_rol=.1, ins_seq_comp=False, check_inv=False, genotype_eval=False,
                out_name = '', fasta_path = None, fasta_id = None, normalize_baseline = False,
                normalize_calls = False):
    """ Run a overlap-based evaluation using the sveval R package (https://github.com/jmonlong/sveval)"""

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

    normalize = normalize_baseline or normalize_calls
    
    # and the fasta
    if fasta_id:
        fasta_name = os.path.basename(fasta_path)
        job.fileStore.readGlobalFile(fasta_id, os.path.join(work_dir, fasta_name))
        # bcftools won't left-align indels over softmasked bases: make sure upper case        
        fa_upper_cmd = ['awk',  'BEGIN{FS=\" \"}{if(!/>/){print toupper($0)}else{print $1}}']
        if fasta_name.endswith('.gz'):
            cmd = [['bgzip', '-d', '-c', fasta_name]]
            if normalize:
                cmd.append(fa_upper_cmd)
            fasta_name = fasta_name[:-3]
            with open(os.path.join(work_dir, fasta_name), 'w') as fasta_file:
                context.runner.call(job, cmd, work_dir = work_dir, outfile=fasta_file)
        elif normalize:
            upper_fasta_name = os.path.splitext(fasta_name)[0] + '_upper.fa'
            with open(os.path.join(work_dir, upper_fasta_name), 'w') as fasta_file:
                context.runner.call(job, fa_upper_cmd + [fasta_name], work_dir = work_dir, outfile=fasta_file)
            fasta_name = upper_fasta_name

    if out_name and not out_name.endswith('_'):
        out_name = '{}_'.format(out_name)
            
    # optional normalization of both calls and truth with bcftools
    if normalize:
        norm_call_vcf_name = '{}calls-norm.vcf.gz'.format(out_name)
        norm_vcfeval_baseline_name = '{}truth-norm.vcf.gz'.format(out_name)
        norm_inputs = []
        if normalize_baseline:
            norm_inputs.append((vcfeval_baseline_name, norm_vcfeval_baseline_name))
        if normalize_calls:
            norm_inputs.append((call_vcf_name, norm_call_vcf_name))
        for vcf_name, norm_name in norm_inputs:
            with open(os.path.join(work_dir, norm_name), 'w') as norm_file:
                # haploid variants throw off bcftools norm --multiallelic +both (TODO: stop making them in vg call)
                norm_cmd = [['bcftools', 'view', vcf_name, '--exclude', 'GT="0" || GT="." || GT="1"']]

                # variants need to be broken up with -both to insure they're fully left-aligned
                norm_cmd.append(['bcftools', 'norm', '-',  '--fasta-ref', fasta_name, '--multiallelic', '-both'])

                # merge up variants back up at the same position
                # (TODO: sveval will need to properly support these cases.  it may end
                #        up being simpler to keep/put them on separate lines at that point)
                norm_cmd.append(['bcftools', 'norm', '-', '--fasta-ref', fasta_name, '--multiallelic', '+both',
                                 '--output-type' ,'z'])
                context.runner.call(job, norm_cmd, work_dir = work_dir, outfile=norm_file)
                
            # bcftools norm --multiallelic -both can apparently unsort the vcf, so we sort it before indexing...
            sort_vcf(job, context.runner, os.path.join(work_dir, norm_name),
                     os.path.join(work_dir, norm_name + '_sorted.vcf'))
            shutil.move(os.path.join(work_dir, norm_name + '_sorted.vcf'), os.path.join(work_dir, norm_name[:-3]))
            context.runner.call(job, ['bgzip', '--force', norm_name[:-3]], work_dir = work_dir)
            context.runner.call(job, ['tabix', '--force', '--preset', 'vcf', norm_name], work_dir = work_dir)
        if normalize_calls:
            call_vcf_name = norm_call_vcf_name
        if normalize_baseline:
            vcfeval_baseline_name = norm_vcfeval_baseline_name

    ## Run sveval R package
    summary_name = '{}sv_accuracy.tsv'.format(out_name)
    sveval_cmd = 'sveval::svevalOl("{}", "{}"'.format(call_vcf_name,
                                                      vcfeval_baseline_name)
    sveval_cmd += ', outfile="{}"'.format(summary_name)
    sveval_cmd += ', min.cov={}'.format(sv_overlap)
    sveval_cmd += ', out.bed.prefix="{}"'.format(out_name)
    if min_sv_len > 0:
        sveval_cmd += ', min.size={}'.format(min_sv_len)
    if max_sv_len < sys.maxint:
        sveval_cmd += ', max.size={}'.format(max_sv_len)
    sveval_cmd += ', max.ins.dist={}'.format(ins_ref_len)
    sveval_cmd += ', min.del.rol={}'.format(del_min_rol)
    if bed_id:
        sveval_cmd += ', bed.regions="{}"'.format(regions_bed_name)
        sveval_cmd += ', bed.regions.ol={}'.format(sv_region_overlap)
    if ins_seq_comp:
        sveval_cmd += ', ins.seq.comp=TRUE'
    if check_inv:
        sveval_cmd += ', check.inv=TRUE'
    if genotype_eval:
        sveval_cmd += ', geno.eval=TRUE, stitch.hets=TRUE, merge.hets=TRUE'
    sveval_cmd += ')'
    r_cmd_file = 'sveval.R'
    with open(os.path.join(work_dir, r_cmd_file), 'w') as r_file:
        r_file.write(sveval_cmd + '\n')
    context.runner.call(job, ['R', '-f', r_cmd_file], work_dir=work_dir)
    summary_id = context.write_output_file(job, os.path.join(work_dir, summary_name))
    
    # tar up some relevant data
    tar_dir = os.path.join(work_dir, '{}sv_evaluation'.format(out_name))
    os.makedirs(tar_dir)
    for dir_file in os.listdir(work_dir):
        if any(dir_file.endswith(ext) for ext in ['.bed', '.tsv', '.vcf.gz', '.vcf.gz.tbi', '.pdf']):
            shutil.copy2(os.path.join(work_dir, dir_file), os.path.join(tar_dir, dir_file))
    context.runner.call(job, ['tar', 'czf', os.path.basename(tar_dir) + '.tar.gz', os.path.basename(tar_dir)],
                        work_dir = work_dir)
    archive_id = context.write_output_file(job, os.path.join(work_dir, tar_dir + '.tar.gz'))

    # Read and return total F1 etc (used in vg_calleval.py)
    results = {}
    with open(os.path.join(work_dir, summary_name)) as summary_file:
        summary_lines = [line for line in summary_file]
    if len(summary_lines) == 1:
        # if no results, we assume there were no SVs and just fill with 0s
        summary_lines.append('\t'.join(['0'] * len(summary_lines[0])))
    headers = summary_lines[0].rstrip().split('\t')
    total_res = summary_lines[1].rstrip().split('\t')
    for idx in range(len(headers)):
        if headers[idx] in ['precision', 'recall', 'F1']:
            try:
                results[headers[idx]] = total_res[idx]
                results[headers[idx]] = float(total_res[idx])
            except:
                results[headers[idx]] = 'error'
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
            
            importer = AsyncImporter(toil)
            
            # Upload local files to the remote IO Store
            vcfeval_baseline_id = importer.load(options.vcfeval_baseline)
            call_vcf_id = importer.load(options.call_vcf)
            vcfeval_baseline_tbi_id = importer.load(options.vcfeval_baseline + '.tbi', wait_on = vcfeval_baseline_id)
            call_tbi_id = importer.load(options.call_vcf + '.tbi', wait_on = call_vcf_id)
            fasta_id = importer.load(options.vcfeval_fasta) if options.vcfeval_fasta else None
            bed_id = importer.load(options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None

            importer.wait()
            vcfeval_baseline_id = importer.resolve(vcfeval_baseline_id)
            call_vcf_id = importer.resolve(call_vcf_id)
            vcfeval_baseline_tbi_id = importer.resolve(vcfeval_baseline_tbi_id)
            call_tbi_id = importer.resolve(call_tbi_id)
            fasta_id = importer.resolve(fasta_id)
            bed_id = importer.resolve(bed_id)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)

            # extract the sample
            if options.vcfeval_sample:
                call_sample_job = init_job.addChildJobFn(run_make_control_vcfs, context, call_vcf_id,
                                                         os.path.basename(options.call_vcf),
                                                         call_tbi_id, options.vcfeval_sample, pos_only = True,
                                                         no_filter_if_sample_not_found = True,
                                                         cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                                         disk=context.config.vcfeval_disk)
                call_vcf_id = call_sample_job.rv(0)
                call_tbi_id = call_sample_job.rv(1)
                
                truth_sample_job = init_job.addChildJobFn(run_make_control_vcfs, context, vcfeval_baseline_id,
                                                          os.path.basename(options.vcfeval_baseline),
                                                          vcfeval_baseline_tbi_id, options.vcfeval_sample, pos_only = True,
                                                          no_filter_if_sample_not_found = True,
                                                          cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                                          disk=context.config.vcfeval_disk)
                vcfeval_baseline_id = truth_sample_job.rv(0)
                vcfeval_baseline_tbi_id = truth_sample_job.rv(1)
            
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
                                       options.min_sv_len, options.max_sv_len,
                                       options.sv_overlap, options.sv_region_overlap,
                                       bed_id,
                                       ins_ref_len=options.ins_max_gap,
                                       del_min_rol=options.del_min_rol,
                                       ins_seq_comp=options.ins_seq_comp,
                                       check_inv=options.check_inv,
                                       genotype_eval=options.genotype_eval,
                                       fasta_path=options.vcfeval_fasta,
                                       fasta_id=fasta_id,
                                       normalize_baseline=options.normalize or options.normalize_baseline,
                                       normalize_calls=options.normalize or options.normalize_calls,
                                       cores=context.config.vcfeval_cores, memory=context.config.vcfeval_mem,
                                       disk=context.config.vcfeval_disk)
                init_job.addFollowOn(sv_job)

            # Run the job
            toil.start(init_job)
        else:
            toil.restart()

    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
