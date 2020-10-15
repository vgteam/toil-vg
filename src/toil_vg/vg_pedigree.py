#!/usr/bin/env python
"""
vg_pedigree.py: pedigree map and calling pipeline to produce parental-enhanced mapping and calling output.

"""
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, struct, socket, threading
import string
import getpass
import gzip
import pipes
import shlex

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

from toil_vg.vg_common import *
from toil_vg.vg_call import run_concat_vcfs
from toil_vg.vg_map import *
from toil_vg.vg_surject import *
from toil_vg.vg_config import *
from toil_vg.vg_construct import *
from toil_vg.vg_index import index_parse_args
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.pedigree_analysis import *

logger = logging.getLogger(__name__)

def pedigree_subparser(parser):
    """
    Create a subparser for pedigree workflow.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("proband_name", type=str,
                        help="sample name of proband or sample of interest (ex HG002)")
    parser.add_argument("maternal_name", type=str,
                        help="sample name of mother to the proband (ex HG004)")
    parser.add_argument("paternal_name", type=str,
                        help="sample name of father to the proband (ex HG003)")
    parser.add_argument("--sibling_names", nargs='+', type=str, default=None,
                        help="sample names of siblings to the proband. Optional.")
    parser.add_argument("--kmer_size", type=int,
                        help="size of kmers to use in gcsa-kmer mapping mode")
    parser.add_argument("--ref_fasta", type=make_url, default=None,
                        help="Path to file with reference fasta.")
    parser.add_argument("--ref_fasta_index", type=make_url, default=None,
                        help="Path to file with reference fasta index.")
    parser.add_argument("--ref_fasta_dict", type=make_url, default=None,
                        help="Path to file with reference fasta dict index.")
    parser.add_argument("--path_list", type=make_url, default=None,
                        help="Path to file with path names of the graph indexes used in the pedigree workflow. One path name per line.")
    parser.add_argument("--ped_file", type=make_url, default=None,
                        help="Path to file containing pedigree definition in .ped format.")
    parser.add_argument("--id_ranges", type=make_url, default=None,
                        help="Path to file with node id ranges for each chromosome in BED format.")
    parser.add_argument("--genetic_map", type=make_url, default=None,
                        help="Path to .tar file containing genetic crossover map files for whatshap phasing.")
    parser.add_argument("--indel_realign_bams", action="store_true", default=False,
                        help="run gatk indel realign on final cohort bams.")
    parser.add_argument("--snpeff_annotation", action="store_true", default=False,
                        help="run snpeff annotation on the final cohort vcf.")
    parser.add_argument("--snpeff_database", type=make_url, default=None,
                        help="Path to .gz file containing snpeff database files for snpeff annotation.")
    parser.add_argument("--run_dragen", action="store_true", default=False,
                        help="run the Illumina Dragen module to do variant calling (NIH Biowulf only).")
    parser.add_argument("--dragen_ref_index_name", type=str, default=None,
                        help="basename of a dragen reference index directory (ex hs37d5_v7) (NIH Biowulf only)")
    parser.add_argument("--udp_data_dir", type=str, default=None,
                        help="basename of a udp data directory which is used to ferry input and output files to and from the dragen module (NIH Biowulf only).")
    parser.add_argument("--helix_username", type=str, default=None,
                        help="username that's used for accessing the Dragen module from NIH's Helix server (NIH Biowulf only)")
    parser.add_argument("--fastq_proband", nargs='+', type=make_url,
                        help="Proband input fastq(s) (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fastq_maternal", nargs='+', type=make_url,
                        help="Maternal input fastq(s) (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fastq_paternal", nargs='+', type=make_url,
                        help="Paternal input fastq(s) (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fastq_siblings", nargs='+', type=make_url, default=None,
                        help="Sibling input fastq(s) (possibly compressed), two are allowed, one for each mate per sibling.\
                            Sibling read-pairs must be input adjacent to eachother. Must follow same order as input to\
                            --sibling_names argument.")
    parser.add_argument("--gam_input_reads_proband", type=make_url, default=None,
                        help="Input reads of proband in GAM format")
    parser.add_argument("--gam_input_reads_maternal", type=make_url, default=None,
                        help="Input reads of mother in GAM format")
    parser.add_argument("--gam_input_reads_paternal", type=make_url, default=None,
                        help="Input reads of father in GAM format")
    parser.add_argument("--gam_input_reads_siblings", nargs='+', type=make_url, default=None,
                        help="Input reads of sibling(s) in GAM format. Must follow same order as input to\
                            --sibling_names argument.")
    parser.add_argument("--bam_input_reads_proband", type=make_url, default=None,
                        help="Input reads of proband in BAM format")
    parser.add_argument("--bam_input_reads_maternal", type=make_url, default=None,
                        help="Input reads of mother in BAM format")
    parser.add_argument("--bam_input_reads_paternal", type=make_url, default=None,
                        help="Input reads of father in BAM format")
    parser.add_argument("--bam_input_reads_siblings", nargs='+', type=make_url, default=None,
                        help="Input reads of sibling(s) in BAM format. Must follow same order as input to\
                            --sibling_names argument.")
    

    # Add common indexing options shared with vg_index
    index_parse_args(parser)

    # Add mapping index options
    map_parse_index_args(parser)

    # Add pedigree options shared only with map
    pedigree_parse_args(parser)
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)
    
    # Add common docker options
    add_container_tool_parse_args(parser)
    
    # Add common analysis options
    pedigree_analysis_parse_args(parser)

def pedigree_parse_args(parser, stand_alone = False):
    """
    Define pedigree arguments shared with map
    """

    parser.add_argument("--fq_split_cores", type=int,
                        help="number of threads used to split input FASTQs")
    parser.add_argument("--single_reads_chunk", action="store_true", default=False,
                        help="do not split reads into chunks")
    parser.add_argument("--reads_per_chunk", type=int,
                        help="number of reads for each mapping job")
    parser.add_argument("--alignment_cores", type=int,
                        help="number of threads during the alignment step")
    parser.add_argument("--interleaved", action="store_true", default=False,
                        help="treat fastq as interleaved read pairs. overrides *_opts")
    parser.add_argument("--map_opts", type=str,
                        help="arguments for vg map (wrapped in \"\")")
    parser.add_argument("--mpmap_opts", type=str,
                        help="arguments for vg mpmap (wrapped in \"\")")
    parser.add_argument("--gaffe_opts", type=str,
                        help="arguments for vg gaffe (wrapped in \"\")")
    parser.add_argument("--bam_output", action="store_true",
                        help="write BAM output directly")
    parser.add_argument("--surject", action="store_true",
                        help="surject output, producing BAM in addition to GAM alignments")
    parser.add_argument("--validate", action="store_true",
                        help="run vg validate on ouput GAMs")
    parser.add_argument("--use_decoys", action="store_true",
                        help="include decoy contigs during parental graph construction")

def pedigree_analysis_parse_args(parser, stand_alone = False):
    """
    Define pedigree arguments related to pedigree analysis
    """
    parser.add_argument("--run_analysis", action="store_true",
                        help="additionally run the pedigree candidate analysis workflow after the vg pedigree workflow completes.")
    # VCFtoShebang options
    parser.add_argument("--bypass", action="store_true", default=False,
        help="Parameter for vcftoshebang.")
    parser.add_argument("--cadd_lines", type=int, default=4985,
        help="Parameter for vcftoshebang.")
    parser.add_argument("--chrom_dir", type=str,
        help="Path to chromosome annotation directory used by vcftoshebang")
    parser.add_argument("--edit_dir", type=str,
        help="Path to directory containing master edit files used by vcftoshebang")
    # CADD options
    parser.add_argument("--split_lines", type=int, default=30000,
        help="Number of lines to chunk the input VCF for CADD processing.")
    parser.add_argument("--genome_build", type=str, default="GRCh37",
        help="Genome annotation version for the CADD engine")
    parser.add_argument("--cadd_data", type=str,
        help="Path to cadd engine data directory")
    # BMTB options
    parser.add_argument("--sibling_genders", nargs='+', type=int,
        help="Gender of each sibling sample. 0 = male, 1 = female. Same sample order as --siblings_bam.")
    parser.add_argument("--sibling_affected", nargs='+', type=int,
        help="Affected status of each sibling sample. 0 = unaffected, 1 = affected. Same sample order as --siblings_bam.")

def validate_pedigree_options(context, options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(options.xg_index is not None, 'All mappers require --xg_index')
    
    if options.mapper == 'map' or options.mapper == 'mpmap':
        require(options.gcsa_index, '--gcsa_index is required for map and mpmap')
    
    if options.mapper == 'gaffe':
        require(options.minimizer_index, '--minimizer_index is required for gaffe')
        require(options.distance_index, '--distance_index is required for gaffe')
        require(options.gbwt_index, '--gbwt_index is required for gaffe')
        require(not options.bam_input_reads, '--bam_input_reads is not supported with gaffe')
        require(not options.interleaved, '--interleaved is not supported with gaffe')
        require(options.fastq is None or len(options.fastq) < 2, 'Multiple --fastq files are not supported with gaffe')
    
    
    require(options.fastq_proband is None or len(options.fastq_proband) in [1, 2], 'Exacty 1 or 2'\
            ' files must be passed with --fastq_proband')
    require(options.fastq_maternal is None or len(options.fastq_maternal) in [1, 2], 'Exacty 1 or 2'\
            ' files must be passed with --fastq_maternal')
    require(options.fastq_paternal is None or len(options.fastq_paternal) in [1, 2], 'Exacty 1 or 2'\
            ' files must be passed with --fastq_paternal')
    require(options.fastq_siblings is None or len(options.fastq_siblings) in [1*len(options.sibling_names), 2*len(options.sibling_names)], 'Exacty 1 or 2 files must be passed per sibling with --fastq_siblings')
    require(options.interleaved == False
            or (options.fastq_proband is None and options.fastq_maternal is None and options.fastq_paternal is None and options.fastq_siblings is None)
            or (len(options.fastq_proband) == 1 and len(options.fastq_maternal) == 1 and len(options.fastq_paternal) == 1 and len(options.fastq_siblings) == len(options.sibling_names)),
            '--interleaved cannot be used when > 1 fastq given for any individual in the pedigree')
    require(sum(map(lambda x : 1 if x else 0, [options.fastq_proband, options.gam_input_reads_proband, options.bam_input_reads_proband])) == 1,
            'reads must be speficied with either --fastq_proband or --gam_input_reads_proband or --bam_input_reads_proband')
    require(sum(map(lambda x : 1 if x else 0, [options.fastq_maternal, options.gam_input_reads_maternal, options.bam_input_reads_maternal])) == 1,
            'reads must be speficied with either --fastq_maternal or --gam_input_reads_maternal or --bam_input_reads_maternal')
    require(sum(map(lambda x : 1 if x else 0, [options.fastq_paternal, options.gam_input_reads_paternal, options.bam_input_reads_paternal])) == 1,
            'reads must be speficied with either --fastq_paternal or --gam_input_reads_paternal or --bam_input_reads_paternal')
    require(options.mapper == 'mpmap' or options.snarls_index is None,
            '--snarls_index can only be used with --mapper mpmap') 
    if options.mapper == 'mpmap':
        require(('-F' in context.config.mpmap_opts or '--output-fmt' in context.config.mpmap_opts) and 'GAM' in context.config.mpmap_opts,
                '-F GAM must be used with mpmap mapper to produce GAM output')
        require(not options.bam_output,
                '--bam_output not currently supported with mpmap mapper')
    require (not options.bam_output or not options.surject,
             '--bam_output cannot be used in combination with --surject')
    require (not options.id_ranges or not options.surject,
             '--surject not currently supported with --id_ranges')
    require (options.snpeff_annotation or options.snpeff_database is None,
             '--snpeff_annotation must be accompanied with --snpeff_database')
    if options.run_dragen:
        require (options.dragen_ref_index_name and options.udp_data_dir and options.helix_username,
             '--run_dragen must be accompanied with --dragen_ref_index_name, --udp_data_dir, and --helix_username {},{},{},{}'.format(options.run_dragen,options.dragen_ref_index_name,options.udp_data_dir,options.helix_username))
    # Requirements for analysis workflow
    if options.run_analysis:
        require(options.chrom_dir, '--chrom_dir is required for analysis workflow')
        require(options.edit_dir, '--edit_dir is required for analysis workflow')
        require(options.cadd_data, '--cadd_data is required for analysis workflow')
        require(len(options.sibling_genders) >= 1, '--sibling_genders needs at least one value for the proband')
        require(len(options.sibling_affected) >= 1, '--sibling_affected needs at least one value for the proband')
    
# Decorator for python process timed retries (https://realpython.com/python-sleep/#adding-a-python-sleep-call-with-decorators)
def sleep(timeout, retry=3):
    def the_real_decorator(function):
        def wrapper(*args, **kwargs):
            retries = 0
            while retries < retry:
                try:
                    value = function(*args, **kwargs)
                    if value is None:
                        return
                except:
                    print(f'Sleeping for {timeout} seconds')
                    time.sleep(timeout)
                    retries += 1
        return wrapper
    return the_real_decorator

def run_gatk_haplotypecaller_gvcf(job, context, sample_name, chr_bam_id, ref_fasta_id,
                                    ref_fasta_index_id, ref_fasta_dict_id, pcr_indel_model="CONSERVATIVE"):

    RealtimeLogger.info("Starting gatk haplotypecalling gvcfs")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample bam for variant calling
    bam_name = os.path.basename(chr_bam_id)
    bam_path = os.path.join(work_dir, bam_name)
    job.fileStore.readGlobalFile(chr_bam_id, bam_path)
    bam_name = os.path.splitext(bam_name)[0]
    
    ref_fasta_name = os.path.basename(ref_fasta_id)
    ref_fasta_path = os.path.join(work_dir, ref_fasta_name)
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    ref_fasta_index_path = os.path.join(work_dir, '{}.fai'.format(ref_fasta_name))
    job.fileStore.readGlobalFile(ref_fasta_index_id, ref_fasta_index_path)
    
    ref_fasta_dict_path = os.path.join(work_dir, '{}.dict'.format(os.path.splitext(ref_fasta_name)[0]))
    job.fileStore.readGlobalFile(ref_fasta_dict_id, ref_fasta_dict_path)
    
    # Extract contig name
    contig_name = re.search('bam_(2[0-2]|1\d|\d|X|Y|MT)', bam_name).group(1)
    
    # Run variant calling commands
    cmd_list = []
    cmd_list.append(['samtools', 'sort', '--threads', str(job.cores), '-n', '-O', 'BAM', os.path.basename(bam_path)])
    cmd_list.append(['samtools', 'fixmate', '-O', 'BAM', '-', '-'])
    cmd_list.append(['samtools', 'sort', '--threads', str(job.cores), '-O', 'BAM', '-'])
    cmd_list.append(['samtools', 'addreplacerg', '-O', 'BAM',
                        '-r', 'ID:1', 
                        '-r', 'LB:lib1',
                        '-r', 'SM:{}'.format(sample_name),
                        '-r', 'PL:illumina',
                        '-r', 'PU:unit1',
                        '-'])
    cmd_list.append(['samtools', 'view', '-@', str(job.cores), '-h', '-O', 'SAM', '-'])
    cmd_list.append(['samtools', 'view', '-@', str(job.cores), '-h', '-O', 'BAM', '-'])
    cmd_list.append(['samtools', 'calmd', '-b', '-', os.path.basename(ref_fasta_path)])
    with open(os.path.join(work_dir, '{}_positionsorted.mdtag.bam'.format(sample_name)), 'wb') as output_samtools_bam:
        context.runner.call(job, cmd_list, work_dir = work_dir, tool_name='samtools', outfile=output_samtools_bam)
    command = ['samtools', 'index', '{}_positionsorted.mdtag.bam'.format(sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='samtools')
    command = ['java', '-Xmx{}'.format(job.memory), '-XX:ParallelGCThreads={}'.format(job.cores), '-jar', '/usr/picard/picard.jar', 'MarkDuplicates',
                'PROGRAM_RECORD_ID=null', 'VALIDATION_STRINGENCY=LENIENT', 'I={}_positionsorted.mdtag.bam'.format(sample_name),
                'O={}.mdtag.dupmarked.bam'.format(sample_name), 'M=marked_dup_metrics.txt']
    with open(os.path.join(work_dir, 'mark_dup_stderr.txt'), 'wb') as outerr_markdupes:
        context.runner.call(job, command, work_dir = work_dir, tool_name='picard', errfile=outerr_markdupes)
    command = ['java', '-Xmx{}'.format(job.memory), '-XX:ParallelGCThreads={}'.format(job.cores), '-jar', '/usr/picard/picard.jar', 'ReorderSam',
                'VALIDATION_STRINGENCY=LENIENT', 'REFERENCE_SEQUENCE={}'.format(os.path.basename(ref_fasta_path)), 'SEQUENCE_DICTIONARY={}'.format(os.path.basename(ref_fasta_dict_path)),
                'INPUT={}.mdtag.dupmarked.bam'.format(sample_name), 'OUTPUT={}_{}.mdtag.dupmarked.reordered.bam'.format(bam_name, sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='picard')
    command = ['samtools', 'index', '{}_{}.mdtag.dupmarked.reordered.bam'.format(bam_name, sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='samtools')
    command = ['gatk', 'HaplotypeCaller',
                '--native-pair-hmm-threads', job.cores,
                '-ERC', 'GVCF',
                '-L', contig_name,
                '--pcr-indel-model', pcr_indel_model,
                '--reference', os.path.basename(ref_fasta_path),
                '--input', '{}_{}.mdtag.dupmarked.reordered.bam'.format(bam_name, sample_name),
                '--output', '{}.{}.rawLikelihoods.gvcf'.format(bam_name, sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='gatk')
    context.runner.call(job, ['bgzip', '{}.{}.rawLikelihoods.gvcf'.format(bam_name, sample_name)],
                        work_dir = work_dir, tool_name='vg')
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', '{}.{}.rawLikelihoods.gvcf.gz'.format(bam_name, sample_name)], work_dir=work_dir)
    
    # Write output to intermediate store
    out_file = os.path.join(work_dir, '{}.{}.rawLikelihoods.gvcf.gz'.format(bam_name, sample_name))
    vcf_file_id = context.write_intermediate_file(job, out_file)
    vcf_index_file_id = context.write_intermediate_file(job, out_file + '.tbi')
    out_bam_file = os.path.join(work_dir, '{}_{}.mdtag.dupmarked.reordered.bam'.format(bam_name, sample_name))
    processed_bam_file_id = context.write_intermediate_file(job, out_bam_file)
    
    return (vcf_file_id, vcf_index_file_id, processed_bam_file_id)

def run_process_chr_bam(job, context, sample_name, chr_bam_id, ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id):
    
    RealtimeLogger.info("Starting bam processing for GATK and Dragen compatibility")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample bam for processing
    bam_name = os.path.basename(chr_bam_id)
    bam_path = os.path.join(work_dir, bam_name)
    job.fileStore.readGlobalFile(chr_bam_id, bam_path)
    bam_name = os.path.splitext(bam_name)[0]
    
    ref_fasta_name = os.path.basename(ref_fasta_id)
    ref_fasta_path = os.path.join(work_dir, ref_fasta_name)
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    ref_fasta_index_path = os.path.join(work_dir, '{}.fai'.format(ref_fasta_name))
    job.fileStore.readGlobalFile(ref_fasta_index_id, ref_fasta_index_path)
    
    ref_fasta_dict_path = os.path.join(work_dir, '{}.dict'.format(os.path.splitext(ref_fasta_name)[0]))
    job.fileStore.readGlobalFile(ref_fasta_dict_id, ref_fasta_dict_path)
    
    # Run variant calling commands
    cmd_list = []
    cmd_list.append(['samtools', 'sort', '--threads', str(job.cores), '-n', '-O', 'BAM', os.path.basename(bam_path)])
    cmd_list.append(['samtools', 'fixmate', '-O', 'BAM', '-', '-'])
    cmd_list.append(['samtools', 'sort', '--threads', str(job.cores), '-O', 'BAM', '-'])
    cmd_list.append(['samtools', 'addreplacerg', '-O', 'BAM',
                        '-r', 'ID:1', 
                        '-r', 'LB:lib1',
                        '-r', 'SM:{}'.format(sample_name),
                        '-r', 'PL:illumina',
                        '-r', 'PU:unit1',
                        '-'])
    cmd_list.append(['samtools', 'view', '-@', str(job.cores), '-h', '-O', 'SAM', '-'])
    cmd_list.append(['samtools', 'view', '-@', str(job.cores), '-h', '-O', 'BAM', '-'])
    cmd_list.append(['samtools', 'calmd', '-b', '-', os.path.basename(ref_fasta_path)])
    with open(os.path.join(work_dir, '{}_positionsorted.mdtag.bam'.format(sample_name)), 'wb') as output_samtools_bam:
        context.runner.call(job, cmd_list, work_dir = work_dir, tool_name='samtools', outfile=output_samtools_bam)
    command = ['samtools', 'index', '{}_positionsorted.mdtag.bam'.format(sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='samtools')
    command = ['java', '-Xmx{}g'.format(int(float(job.memory)/2000000000)), '-XX:ParallelGCThreads={}'.format(job.cores), '-jar', '/usr/picard/picard.jar', 'MarkDuplicates',
                'PROGRAM_RECORD_ID=null', 'VALIDATION_STRINGENCY=LENIENT', 'I={}_positionsorted.mdtag.bam'.format(sample_name),
                'O={}.mdtag.dupmarked.bam'.format(sample_name), 'M=marked_dup_metrics.txt']
    with open(os.path.join(work_dir, 'mark_dup_stderr.txt'), 'wb') as outerr_markdupes:
        context.runner.call(job, command, work_dir = work_dir, tool_name='picard', errfile=outerr_markdupes)
    command = ['java', '-Xmx{}g'.format(int(float(job.memory)/1000000000)), '-XX:ParallelGCThreads={}'.format(job.cores), '-jar', '/usr/picard/picard.jar', 'ReorderSam',
                'VALIDATION_STRINGENCY=LENIENT', 'REFERENCE_SEQUENCE={}'.format(os.path.basename(ref_fasta_path)), 'SEQUENCE_DICTIONARY={}'.format(os.path.basename(ref_fasta_dict_path)),
                'INPUT={}.mdtag.dupmarked.bam'.format(sample_name), 'OUTPUT={}_{}.mdtag.dupmarked.reordered.bam'.format(bam_name, sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='picard')
    
    # Write output to intermediate store
    out_bam_file = os.path.join(work_dir, '{}_{}.mdtag.dupmarked.reordered.bam'.format(bam_name, sample_name))
    processed_bam_file_id = context.write_intermediate_file(job, out_bam_file)
    
    # Delete input files
    job.fileStore.deleteGlobalFile(chr_bam_id)
    
    return processed_bam_file_id

#@sleep(30, retry=20)
@sleep(900, retry=20)
def run_dragen_commands(job, context, command, work_dir):
    """ 
    Helper function for running the Dragen gvcf caller asynchronously
    """
    try:
        context.runner.call(job, command, work_dir = work_dir)
    except:
        raise

def run_dragen_gvcf(job, context, sample_name, merge_bam_id, dragen_ref_index_name, udp_data_dir, helix_username, write_to_outstore=False):
    
    RealtimeLogger.info("Starting Dragen GVCF caller")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample bam for processing
    bam_name = os.path.basename(merge_bam_id)
    bam_path = os.path.join(work_dir, bam_name)
    job.fileStore.readGlobalFile(merge_bam_id, bam_path)
    
    udp_data_dir_path = '{}/usr/{}'.format(udp_data_dir, helix_username)
    dragen_work_dir_path = '/staging/{}/{}'.format(helix_username, sample_name)
    tmp_dir_path = '/staging/{}/tmp'.format(helix_username)
    udp_data_bam_path = '/data/{}/{}_surjected_bams/'.format(udp_data_dir_path, sample_name)
    
    # Make sure directory paths are valid
    assert ' ' not in udp_data_dir_path
    assert ' ' not in dragen_work_dir_path
    assert ' ' not in tmp_dir_path
    assert ' ' not in udp_data_bam_path
    cmd_list = []
    cmd_list.append(['mkdir', '-p', udp_data_bam_path])
    cmd_list.append(['cp', str(bam_path), udp_data_bam_path])
    cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"mkdir -p {}\"'.format(dragen_work_dir_path)])
    cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"mkdir -p {}\"'.format(tmp_dir_path)])
    cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username),
                              '\"' +
                              'dragen -f -r /staging/{}'.format(dragen_ref_index_name) +
                              ' -b /staging/helix/{}/{}_surjected_bams/{}'.format(udp_data_dir_path, sample_name, bam_name) +
                              ' --verbose --bin_memory=50000000000 --enable-map-align false --enable-variant-caller true' +
                              ' --pair-by-name=true --vc-emit-ref-confidence GVCF' +
                              ' --intermediate-results-dir {} --output-directory {} --output-file-prefix {}_dragen_genotyped'.format(tmp_dir_path, dragen_work_dir_path, sample_name) +
                              '\"'])
    cmd_list.append(['mkdir', '/data/{}/{}_dragen_genotyper'.format(udp_data_dir_path, sample_name)])
    cmd_list.append(['chmod', 'ug+rw', '-R', '/data/{}/{}_dragen_genotyper'.format(udp_data_dir_path, sample_name)])
    cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"cp -R {} /staging/helix/{}/{}_dragen_genotyper \"'.format(dragen_work_dir_path, udp_data_dir_path, sample_name)])
    cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"rm -fr {}/\"'.format(dragen_work_dir_path)])
    cmd_list.append(['mv', '/data/{}/{}_dragen_genotyper'.format(udp_data_dir_path, sample_name), '{}_dragen_genotyper'.format(sample_name)])
    cmd_list.append(['rm', '-f', '{}{}'.format(udp_data_bam_path, bam_name)])
    cmd_list.append(['rmdir', '{}'.format(udp_data_bam_path)])
    chain_cmds = [' '.join(p) for p in cmd_list]
    command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' && '.join(chain_cmds))]
    run_dragen_commands(job, context, command, work_dir)
    
    # Write output to intermediate store
    out_gvcf_file = os.path.join(work_dir, '{}_dragen_genotyper/{}/{}_dragen_genotyped.hard-filtered.gvcf.gz'.format(sample_name, sample_name, sample_name))
    if write_to_outstore:
        processed_gvcf_file_id = context.write_output_file(job, out_gvcf_file)
        processed_gvcf_index_file_id = context.write_output_file(job, out_gvcf_file + '.tbi')
    else:
        processed_gvcf_file_id = context.write_intermediate_file(job, out_gvcf_file)
        processed_gvcf_index_file_id = context.write_intermediate_file(job, out_gvcf_file + '.tbi')
    
    return (processed_gvcf_file_id, processed_gvcf_index_file_id)

def run_merge_bams_ped_workflow(job, context, sample_name, bam_ids, indel_realign_bams_name=False, write_to_outstore = False):
    
    RealtimeLogger.info("Starting samtools merge bams")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download the input
    bam_paths = []
    for id_num,bam_id in enumerate(bam_ids):
        bam_path = os.path.join(work_dir, "{}.{}".format(os.path.basename(bam_id),id_num))
        job.fileStore.readGlobalFile(bam_id, bam_path)
        bam_paths.append(os.path.basename(bam_path))
    
    out_file = os.path.join(work_dir, '{}_merged.bam'.format(sample_name))
    if indel_realign_bams_name:
        out_file = os.path.join(work_dir, '{}_merged.indel_realigned.bam'.format(sample_name))
    
    command = ['samtools', 'merge', '-f', '-p', '-c', '--threads', job.cores,
                    os.path.basename(out_file)] + bam_paths
    context.runner.call(job, command, work_dir = work_dir, tool_name='samtools')
    command = ['samtools', 'index', os.path.basename(out_file)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='samtools')
    
    if write_to_outstore:
        merged_bam_file_id = context.write_output_file(job, out_file)
        merged_bam_index_file_id = context.write_output_file(job, out_file + '.bai')
        # Delete input files
        for bam_id in bam_ids:
            job.fileStore.deleteGlobalFile(bam_id)
    else:
        merged_bam_file_id = context.write_intermediate_file(job, out_file)
        merged_bam_index_file_id = context.write_intermediate_file(job, out_file + '.bai')
    
    
    return (merged_bam_file_id, merged_bam_index_file_id)
    
def run_pipeline_call_gvcfs(job, context, options, sample_name, chr_bam_ids, ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                run_dragen=False, dragen_ref_index_name=None, udp_data_dir=None, helix_username=None):
    """
    Call all the chromosomes and return a merged up vcf/tbi pair
    """
    RealtimeLogger.info("Starting gvcf calling pipeline for sample: {}".format(sample_name))
    vcf_ids = []
    tbi_ids = []
    processed_bam_ids = []
    # Write bams to final outstore only if not also doing indel realignment
    if options.indel_realign_bams:
        write_to_outstore = False
    else:
        write_to_outstore = True
    
    # If running without dragen then
    child_job = Job()
    job.addChild(child_job)
    if not run_dragen:
        for chr_bam_id in chr_bam_ids:
            call_job = child_job.addChildJobFn(run_gatk_haplotypecaller_gvcf, context, sample_name, chr_bam_id, ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                            cores=context.config.alignment_cores, memory=context.config.alignment_mem, disk=context.config.alignment_disk)
            vcf_ids.append(call_job.rv(0))
            tbi_ids.append(call_job.rv(1))
            processed_bam_ids.append(call_job.rv(2))
    else:
        for chr_bam_id in chr_bam_ids:
            process_bam_job = child_job.addChildJobFn(run_process_chr_bam, context, sample_name, chr_bam_id, ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                            cores=context.config.alignment_cores, memory="{}G".format(int(re.findall(r'\d+', context.config.alignment_mem)[0])*2), disk="{}G".format(int(re.findall(r'\d+', context.config.alignment_disk)[0])*4)) 
            processed_bam_ids.append(process_bam_job.rv())
    
    # Run merging of crhomosomal bams
    merge_chr_bams_job = child_job.addFollowOnJobFn(run_merge_bams_ped_workflow, context, sample_name, processed_bam_ids, False, write_to_outstore,
                                                            cores=context.config.alignment_cores, memory=context.config.alignment_mem, disk=context.config.alignment_disk)
    
    # Run gvcf concatenation
    # If using the Illumina Dragen module for the NIH Biowulf system, then instead run that on the processed merged bam file
    output_gvcf_id = None
    output_gvcf_index_id = None
    if not run_dragen:
        concat_job = merge_chr_bams_job.addChildJobFn(run_concat_vcfs, context, sample_name, vcf_ids, tbi_ids, write_to_outstore = True)
        output_gvcf_id = concat_job.rv(0)
        output_gvcf_index_id = concat_job.rv(1)
    else:
        dragen_job = merge_chr_bams_job.addChildJobFn(run_dragen_gvcf, context, sample_name, merge_chr_bams_job.rv(0), dragen_ref_index_name, udp_data_dir, helix_username, write_to_outstore = True)
        output_gvcf_id = dragen_job.rv(0)
        output_gvcf_index_id = dragen_job.rv(1)
    
    if options.indel_realign_bams:
        return (merge_chr_bams_job.rv(0), merge_chr_bams_job.rv(1), output_gvcf_id, output_gvcf_index_id, processed_bam_ids)
    else:
        return (merge_chr_bams_job.rv(0), merge_chr_bams_job.rv(1), output_gvcf_id, output_gvcf_index_id)

def run_joint_genotyper(job, context, sample_name, proband_gvcf_id, proband_gvcf_index_id,
                                    maternal_gvcf_id, maternal_gvcf_index_id,
                                    paternal_gvcf_id, paternal_gvcf_index_id,
                                    sibling_call_gvcf_ids, sibling_call_gvcf_index_ids,
                                    ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                    snpeff_annotation=False, run_dragen=False,
                                    dragen_ref_index_name=None, udp_data_dir=None, helix_username=None):

    RealtimeLogger.info("Starting gatk joint calling gvcfs")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample gvcfs for joint genotyping
    proband_gvcf_path = os.path.join(work_dir, os.path.basename(proband_gvcf_id))
    job.fileStore.readGlobalFile(proband_gvcf_id, proband_gvcf_path)
    proband_gvcf_index_path = os.path.join(work_dir, '{}.tbi'.format(os.path.basename(proband_gvcf_id)))
    job.fileStore.readGlobalFile(proband_gvcf_index_id, proband_gvcf_index_path)
    
    maternal_gvcf_path = os.path.join(work_dir, os.path.basename(maternal_gvcf_id))
    job.fileStore.readGlobalFile(maternal_gvcf_id, maternal_gvcf_path)
    maternal_gvcf_index_path = os.path.join(work_dir, '{}.tbi'.format(os.path.basename(maternal_gvcf_id)))
    job.fileStore.readGlobalFile(maternal_gvcf_index_id, maternal_gvcf_index_path)
    
    paternal_gvcf_path = os.path.join(work_dir, os.path.basename(paternal_gvcf_id))
    job.fileStore.readGlobalFile(paternal_gvcf_id, paternal_gvcf_path)
    paternal_gvcf_index_path = os.path.join(work_dir, '{}.tbi'.format(os.path.basename(paternal_gvcf_id)))
    job.fileStore.readGlobalFile(paternal_gvcf_index_id, paternal_gvcf_index_path)
    
    ref_fasta_name = os.path.basename(ref_fasta_id)
    ref_fasta_path = os.path.join(work_dir, ref_fasta_name)
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    ref_fasta_index_path = os.path.join(work_dir, '{}.fai'.format(ref_fasta_name))
    job.fileStore.readGlobalFile(ref_fasta_index_id, ref_fasta_index_path)
    
    ref_fasta_dict_path = os.path.join(work_dir, '{}.dict'.format(os.path.splitext(ref_fasta_name)[0]))
    job.fileStore.readGlobalFile(ref_fasta_dict_id, ref_fasta_dict_path)
    
    sibling_options_list = []
    if sibling_call_gvcf_ids is not None and sibling_call_gvcf_index_ids is not None:
        for sibling_gvcf_id, sibling_call_gvcf_index in zip(sibling_call_gvcf_ids,sibling_call_gvcf_index_ids):
            sibling_gvcf_path = os.path.join(work_dir, os.path.basename(sibling_gvcf_id))
            job.fileStore.readGlobalFile(sibling_gvcf_id, sibling_gvcf_path)
            sibling_gvcf_index_path = os.path.join(work_dir, '{}.tbi'.format(os.path.basename(sibling_gvcf_id)))
            job.fileStore.readGlobalFile(sibling_call_gvcf_index, sibling_gvcf_index_path)
            if not run_dragen:
                sibling_options_list += ['-V', os.path.basename(sibling_gvcf_path)]
            else:
                sibling_options_list += [' --variant /staging/helix/{}/usr/{}/{}_cohort_gvcfs/{}'.format(udp_data_dir, helix_username, sample_name, os.path.basename(sibling_gvcf_path))]
    
    out_file = None
    if not run_dragen:        
        # Run gatk variant calling commands
        command = ['gatk', 'CombineGVCFs',
                    '--reference', os.path.basename(ref_fasta_path),
                    '-V', os.path.basename(maternal_gvcf_path),
                    '-V', os.path.basename(paternal_gvcf_path),
                    '-V', os.path.basename(proband_gvcf_path)]
        command += sibling_options_list
        command += ['--output', '{}_trio.combined.gvcf'.format(sample_name)]
        context.runner.call(job, command, work_dir = work_dir, tool_name='gatk')
        command = ['gatk', 'GenotypeGVCFs',
                    '--reference', os.path.basename(ref_fasta_path),
                    '--variant', '{}_trio.combined.gvcf'.format(sample_name),
                    '--output', '{}_trio.jointgenotyped.unnormalized.vcf'.format(sample_name)]
        context.runner.call(job, command, work_dir = work_dir, tool_name='gatk')
        command = ['bcftools', 'norm', '-m-both', '--threads', job.cores,
                        '-o', '{}_trio.jointgenotyped.bcftools_normalized.vcf'.format(sample_name), '{}_trio.jointgenotyped.unnormalized.vcf'.format(sample_name)]
        context.runner.call(job, command, work_dir = work_dir, tool_name='bcftools')
        context.runner.call(job, ['gatk', 'SelectVariants',
                                    '-R', os.path.basename(ref_fasta_path),
                                    '--remove-unused-alternates',
                                    '--exclude-non-variants',
                                    '-V', '{}_trio.jointgenotyped.bcftools_normalized.vcf'.format(sample_name),
                                    '-O', '{}_trio.jointgenotyped.vcf'.format(sample_name)],
                            work_dir = work_dir, tool_name='gatk')
        context.runner.call(job, ['bgzip', '{}_trio.jointgenotyped.vcf'.format(sample_name)],
                            work_dir = work_dir, tool_name='vg')
        context.runner.call(job, ['tabix', '-f', '-p', 'vcf', '{}_trio.jointgenotyped.vcf.gz'.format(sample_name)], work_dir=work_dir)
        # checkpoint to out store
        out_file = os.path.join(work_dir, '{}_trio.jointgenotyped.vcf.gz'.format(sample_name))
    else:
        # Run dragen variant calling commands
        udp_data_dir_path = '{}/usr/{}'.format(udp_data_dir, helix_username)
        joint_genotype_dragen_work_dir_path = '/staging/{}/output_cohort_joint_call_{}'.format(helix_username, sample_name)
        tmp_dir_path = '/staging/{}/tmp'.format(helix_username)
        udp_data_gvcf_path = '/data/{}/{}_cohort_gvcfs/'.format(udp_data_dir_path, sample_name)
        
        # Make sure directory paths are valid
        assert ' ' not in udp_data_dir_path
        assert ' ' not in joint_genotype_dragen_work_dir_path
        assert ' ' not in tmp_dir_path
        assert ' ' not in udp_data_gvcf_path
        
        context.runner.call(job, ['mkdir', '-p', udp_data_gvcf_path], work_dir = work_dir)
        command = ['cp', os.path.basename(maternal_gvcf_path), os.path.basename(paternal_gvcf_path), os.path.basename(proband_gvcf_path)]
        if sibling_call_gvcf_ids is not None and sibling_call_gvcf_index_ids is not None:
            for sibling_gvcf_id in sibling_call_gvcf_ids:
                sibling_gvcf_path = os.path.join(work_dir, os.path.basename(sibling_gvcf_id))
                command += [os.path.basename(sibling_gvcf_path)]
        command += [udp_data_gvcf_path]
        context.runner.call(job, command, work_dir = work_dir)
        cmd_list = []
        cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"mkdir -p {}\"'.format(joint_genotype_dragen_work_dir_path)])
        cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"mkdir -p {}\"'.format(tmp_dir_path)])
        cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username),
                                  '\"' +
                                  'dragen -f -r /staging/{}'.format(dragen_ref_index_name) +
                                  ' --enable-joint-genotyping true --intermediate-results-dir {}'.format(tmp_dir_path) +
                                  ' --output-directory {} --output-file-prefix cohort_joint_genotyped_{}'.format(joint_genotype_dragen_work_dir_path, sample_name) +
                                  ' --variant /staging/helix/{}/{}_cohort_gvcfs/{}'.format(udp_data_dir_path, sample_name, os.path.basename(maternal_gvcf_path)) +
                                  ' --variant /staging/helix/{}/{}_cohort_gvcfs/{}'.format(udp_data_dir_path, sample_name, os.path.basename(paternal_gvcf_path)) +
                                  ' --variant /staging/helix/{}/{}_cohort_gvcfs/{}'.format(udp_data_dir_path, sample_name, os.path.basename(proband_gvcf_path)) +
                                  ' '.join(sibling_options_list) +
                                  '\"'])
        cmd_list.append(['mkdir', '/data/{}/{}_dragen_joint_genotyper'.format(udp_data_dir_path, sample_name)])
        cmd_list.append(['chmod', 'ug+rw', '-R', '/data/{}/{}_dragen_joint_genotyper'.format(udp_data_dir_path, sample_name)])
        cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"cp -R {} /staging/helix/{}/{}_dragen_joint_genotyper \"'.format(joint_genotype_dragen_work_dir_path, udp_data_dir_path, sample_name)])
        cmd_list.append(['ssh', '{}@helix.nih.gov'.format(helix_username), 'ssh', '{}@udpdragen01.nhgri.nih.gov'.format(helix_username), '\"rm -fr {}/\"'.format(joint_genotype_dragen_work_dir_path)])
        cmd_list.append(['mv', '/data/{}/{}_dragen_joint_genotyper'.format(udp_data_dir_path, sample_name), '{}_dragen_joint_genotyper'.format(sample_name)])
        cmd_list.append(['rm', '-f', '{}{}'.format(udp_data_gvcf_path, os.path.basename(maternal_gvcf_path))])
        cmd_list.append(['rm', '-f', '{}{}'.format(udp_data_gvcf_path, os.path.basename(paternal_gvcf_path))])
        cmd_list.append(['rm', '-f', '{}{}'.format(udp_data_gvcf_path, os.path.basename(proband_gvcf_path))])
        if sibling_call_gvcf_ids is not None and sibling_call_gvcf_index_ids is not None:
            for sibling_gvcf_id in sibling_call_gvcf_ids:
                sibling_gvcf_path = os.path.join(work_dir, os.path.basename(sibling_gvcf_id))
                cmd_list.append(['rm', '-f', '{}{}'.format(udp_data_gvcf_path, os.path.basename(sibling_gvcf_path))])
        cmd_list.append(['rmdir', '{}'.format(udp_data_gvcf_path)])
        chain_cmds = [' '.join(p) for p in cmd_list]
        command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' && '.join(chain_cmds))]
        run_dragen_commands(job, context, command, work_dir)
        out_file = os.path.join(work_dir, '{}_dragen_joint_genotyper/output_cohort_joint_call_{}/cohort_joint_genotyped_{}.vcf.gz'.format(sample_name, sample_name, sample_name))
    
    if snpeff_annotation:
        joint_vcf_file_id = context.write_output_file(job, out_file)
        joint_vcf_index_file_id = context.write_output_file(job, out_file + '.tbi')
    else:
        joint_vcf_file_id = context.write_output_file(job, out_file)
        joint_vcf_index_file_id = context.write_output_file(job, out_file + '.tbi')
    
    return (joint_vcf_file_id, joint_vcf_index_file_id)

def run_split_jointcalled_vcf(job, context, joint_called_vcf_id, joint_called_vcf_index_id, proband_name, maternal_name, paternal_name, contigs_list, filter_parents=False):
    
    RealtimeLogger.info("Starting split joint-called trio vcf")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    joint_vcf_name = os.path.basename(joint_called_vcf_id)
    joint_vcf_path = os.path.join(work_dir, joint_vcf_name)
    job.fileStore.readGlobalFile(joint_called_vcf_id, joint_vcf_path)
    
    joint_vcf_index_path = os.path.join(work_dir, '{}.tbi'.format(joint_vcf_name))
    job.fileStore.readGlobalFile(joint_called_vcf_index_id, joint_vcf_index_path)
    
    contig_vcf_pair_ids = {}
    for contig_id in contigs_list:
        # Define temp file for our output
        output_file = os.path.join(work_dir, "{}_trio_{}.vcf.gz".format(proband_name, contig_id))
        
        # Open the file stream for writing
        with open(output_file, "wb") as contig_vcf_file:
            command = []
            if filter_parents == True and contig_id == "MT":
                command = ['bcftools', 'view', '-O', 'z', '-r', contig_id, '-s', maternal_name, joint_vcf_name]
            elif filter_parents == True and contig_id == "Y":
                command = ['bcftools', 'view', '-O', 'z', '-r', contig_id, '-s', paternal_name, joint_vcf_name]
            elif filter_parents == True:
                command = ['bcftools', 'view', '-O', 'z', '-r', contig_id, '-s', '{},{}'.format(maternal_name,paternal_name), joint_vcf_name]
            else:
                command = ['bcftools', 'view', '-O', 'z', '-r', contig_id, joint_vcf_name]
            context.runner.call(job, command, work_dir = work_dir, tool_name='bcftools', outfile=contig_vcf_file)
        contig_vcf_pair_ids[contig_id] = context.write_intermediate_file(job, output_file)
    
    return contig_vcf_pair_ids

def run_whatshap_phasing(job, context, contig_vcf_id, contig_name, proband_name, maternal_name, paternal_name,
                            proband_bam_id, proband_bam_index_id,
                            maternal_bam_id, maternal_bam_index_id,
                            paternal_bam_id, paternal_bam_index_id,
                            ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                            ped_file_id, genetic_map_id):
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    contig_vcf_path =  os.path.join(work_dir, os.path.basename(contig_vcf_id))
    job.fileStore.readGlobalFile(contig_vcf_id, contig_vcf_path)
     
    proband_bam_name = os.path.basename(proband_bam_id)
    proband_bam_path = os.path.join(work_dir, proband_bam_name)
    job.fileStore.readGlobalFile(proband_bam_id, proband_bam_path)
    proband_bam_index_path = os.path.join(work_dir, '{}.bai'.format(proband_bam_name))
    job.fileStore.readGlobalFile(proband_bam_index_id, proband_bam_index_path)
    
    maternal_bam_name = os.path.basename(maternal_bam_id)
    maternal_bam_path = os.path.join(work_dir, maternal_bam_name)
    job.fileStore.readGlobalFile(maternal_bam_id, maternal_bam_path)
    maternal_bam_index_path = os.path.join(work_dir, '{}.bai'.format(maternal_bam_name))
    job.fileStore.readGlobalFile(maternal_bam_index_id, maternal_bam_index_path)
    
    paternal_bam_name = os.path.basename(paternal_bam_id)
    paternal_bam_path = os.path.join(work_dir, paternal_bam_name)
    job.fileStore.readGlobalFile(paternal_bam_id, paternal_bam_path)
    paternal_bam_index_path = os.path.join(work_dir, '{}.bai'.format(paternal_bam_name))
    job.fileStore.readGlobalFile(paternal_bam_index_id, paternal_bam_index_path)
    
    ref_fasta_name = os.path.basename(ref_fasta_id)
    ref_fasta_path = os.path.join(work_dir, ref_fasta_name)
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    ref_fasta_index_path = os.path.join(work_dir, '{}.fai'.format(ref_fasta_name))
    job.fileStore.readGlobalFile(ref_fasta_index_id, ref_fasta_index_path)
    
    ref_fasta_dict_path = os.path.join(work_dir, '{}.dict'.format(os.path.splitext(ref_fasta_name)[0]))
    job.fileStore.readGlobalFile(ref_fasta_dict_id, ref_fasta_dict_path)
    
    ped_file_path = os.path.join(work_dir, os.path.basename(ped_file_id))
    job.fileStore.readGlobalFile(ped_file_id, ped_file_path)
    
    # Run trio-based mendelian correction
    context.runner.call(job, ['gatk', '--java-options', '-Xmx{}g'.format(int(float(job.memory)/1000000000)), 'IndexFeatureFile', '-F', os.path.basename(contig_vcf_path)], work_dir = work_dir, tool_name='gatk')
    context.runner.call(job, ['gatk', '--java-options', '-Xmx{}g'.format(int(float(job.memory)/1000000000)), 'CalculateGenotypePosteriors', 
                              '-V', os.path.basename(contig_vcf_path),
                              '-O', '{}_cohort_{}.gatk_mendelian_corrected.vcf.gz'.format(proband_name, contig_name),
                              '-ped', os.path.basename(ped_file_path), '--skip-population-priors'], work_dir = work_dir, tool_name='gatk')
    # Run eagle phasing
    whatshap_input_file = '{}_cohort_{}.gatk_mendelian_corrected.vcf.gz'.format(proband_name, contig_name)
    if contig_name not in ['Y', 'MT', 'ABOlocus']:
        # Use eagle phasing
        whatshap_input_file = '{}_cohort_{}.gatk_mendelian_corrected.eagle_phased.vcf.gz'.format(proband_name, contig_name)
        context.runner.call(job, ['wget', 'ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz'], work_dir = work_dir)
        with open(os.path.join(work_dir, 'human_g1k_v37.fasta'), "wb") as ref_fasta_file:
            context.runner.call(job, ['gzip', '-d', '-c', 'Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz'], work_dir = work_dir, outfile = ref_fasta_file)
        context.runner.call(job, ['samtools', 'faidx', 'human_g1k_v37.fasta'], work_dir = work_dir, tool_name='samtools')
        if contig_name in ['X']:
            # Run eagle phasing on X chromsome
            context.runner.call(job, ['wget', 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'], work_dir = work_dir)
            context.runner.call(job, ['wget', 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release//20130502/ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz.tbi'], work_dir = work_dir)
            cmd_list = []
            cmd_list.append(['bcftools', 'view', '--no-version', '-Ou', '-c', '2', 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'])
            cmd_list.append(['bcftools', 'norm', '--no-version', '-Ou', '-m', '-any'])
            cmd_list.append(['bcftools', 'norm', '--no-version', '-Ob', '-o', 'ALL.chrX.phase3_integrated.20130502.genotypes.bcf', '-d', 'none', '-f', 'human_g1k_v37.fasta'])
            context.runner.call(job, cmd_list, work_dir = work_dir, tool_name='bcftools')
            context.runner.call(job, ['bcftools', 'index', '-f', 'ALL.chrX.phase3_integrated.20130502.genotypes.bcf'], work_dir = work_dir, tool_name='bcftools')
            context.runner.call(job, ['/usr/src/app/eagle', '--outputUnphased', '--geneticMapFile', '/usr/src/app/genetic_map_hg19_withX.txt.gz', '--outPrefix', '{}_cohort_{}.gatk_mendelian_corrected.eagle_phased'.format(proband_name, contig_name),
                                '--numThreads', job.cores, '--vcfRef', 'ALL.chrX.phase3_integrated.20130502.genotypes.bcf', '--vcfTarget', '{}_cohort_{}.gatk_mendelian_corrected.vcf.gz'.format(proband_name, contig_name),
                                '--chrom', contig_name], work_dir = work_dir, tool_name='eagle')
        else:
            # Run eagle phasing on autosomal chromosome
            context.runner.call(job, ['wget', 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(contig_name)], work_dir = work_dir)
            context.runner.call(job, ['wget', 'ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi'.format(contig_name)], work_dir = work_dir)
            cmd_list = []
            cmd_list.append(['bcftools', 'view', '--no-version', '-Ou', '-c', '2', 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(contig_name)])
            cmd_list.append(['bcftools', 'norm', '--no-version', '-Ou', '-m', '-any'])
            cmd_list.append(['bcftools', 'norm', '--no-version', '-Ob', '-o', 'ALL.chr{}.phase3_integrated.20130502.genotypes.bcf'.format(contig_name), '-d', 'none', '-f', 'human_g1k_v37.fasta'])
            context.runner.call(job, cmd_list, work_dir = work_dir, tool_name='bcftools')
            context.runner.call(job, ['bcftools', 'index', '-f', 'ALL.chr{}.phase3_integrated.20130502.genotypes.bcf'.format(contig_name)], work_dir = work_dir, tool_name='bcftools')
            context.runner.call(job, ['/usr/src/app/eagle', '--outputUnphased', '--geneticMapFile', '/usr/src/app/genetic_map_hg19_withX.txt.gz', '--outPrefix', '{}_cohort_{}.gatk_mendelian_corrected.eagle_phased'.format(proband_name, contig_name),
                                '--numThreads', job.cores, '--vcfRef', 'ALL.chr{}.phase3_integrated.20130502.genotypes.bcf'.format(contig_name), '--vcfTarget', '{}_cohort_{}.gatk_mendelian_corrected.vcf.gz'.format(proband_name, contig_name),
                                '--chrom', contig_name], work_dir = work_dir, tool_name='eagle')
        
    
    # Run whatshap phasing
    command = ['whatshap', 'phase', '--reference', os.path.basename(ref_fasta_path), '--indels', '--ped', os.path.basename(ped_file_path)]
    
    if bool(genetic_map_id) == True:
        genetic_map_path = os.path.join(work_dir, os.path.basename(genetic_map_id))
        job.fileStore.readGlobalFile(genetic_map_id, genetic_map_path)
        context.runner.call(job, ['tar', '-xvf', os.path.basename(genetic_map_path)], work_dir = work_dir, tool_name='whatshap')
        if contig_name not in ['X','Y', 'MT', 'ABOlocus']:
            command += ['--genmap', 'genetic_map_GRCh37/genetic_map_chr{}_combined_b37.txt'.format(contig_name), '--chromosome', contig_name]
        elif contig_name == 'X':
            command += ['--genmap', 'genetic_map_GRCh37/genetic_map_chrX_nonPAR_combined_b37.txt', '--chromosome', 'X']

    command += ['-o', '{}_cohort_{}.phased.vcf'.format(proband_name, contig_name), whatshap_input_file,
                    proband_bam_name, maternal_bam_name, paternal_bam_name]
    context.runner.call(job, command, work_dir = work_dir, tool_name='whatshap')
    
    # Filter for phased parental genotypes ONLY
    context.runner.call(job, ['bgzip', '{}_cohort_{}.phased.vcf'.format(proband_name, contig_name)], work_dir = work_dir, tool_name='whatshap')
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', '{}_cohort_{}.phased.vcf.gz'.format(proband_name, contig_name)], work_dir=work_dir)
    #command = []
    #command.append(['bcftools', 'view', '-Oz', '-r', contig_name, '-s', '{},{}'.format(maternal_name,paternal_name), '{}_cohort_{}.phased.vcf.gz'.format(proband_name, contig_name)])
    #command.append(['bcftools', 'view', '-p', '-Oz', '-'])
    command = ['bcftools', 'view', '-Oz', '-r', contig_name, '-s', '{},{}'.format(maternal_name,paternal_name), '{}_cohort_{}.phased.vcf.gz'.format(proband_name, contig_name)]
    output_file = os.path.join(work_dir, '{}.vcf.gz'.format(contig_name))
    with open(output_file, "wb") as contig_vcf_file:
            context.runner.call(job, command, work_dir = work_dir, tool_name='bcftools', outfile=contig_vcf_file)
    
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', '{}.vcf.gz'.format(contig_name)], work_dir=work_dir)
    # Write output to intermediate store
    out_file = os.path.join(work_dir, '{}.vcf.gz'.format(contig_name))
    phased_vcf_file_id = context.write_output_file(job, out_file)
    phased_vcf_index_file_id = context.write_output_file(job, out_file + '.tbi')
    phased_vcf_file_name = os.path.basename(out_file)
    
    # Delete input files
    job.fileStore.deleteGlobalFile(contig_vcf_id)
    
    return (phased_vcf_file_id, phased_vcf_index_file_id, phased_vcf_file_name)

def run_collect_concat_vcfs(job, context, vcf_file_id, vcf_index_file_id):
    inputVCFFileIDs = []
    inputVCFNames = []
    inputTBIFileIDs = []
    
    inputVCFFileIDs.append(vcf_file_id)
    inputVCFNames.append(os.path.basename(vcf_file_id))
    inputTBIFileIDs.append(vcf_index_file_id)
    return (inputVCFFileIDs, inputVCFNames, inputTBIFileIDs)

def run_pipeline_construct_parental_graphs(job, context, options, joint_called_vcf_id, joint_called_vcf_index_id, proband_name, maternal_name, paternal_name, 
                                            proband_bam_id, proband_bam_index_id,
                                            maternal_bam_id, maternal_bam_index_id,
                                            paternal_bam_id, paternal_bam_index_id,
                                            ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                            path_list_id, ped_file_id, genetic_map_id):
    
    # we make a sub job tree so that all phasing and graph construction is encapsulated in a top-level job
    split_job = Job()
    phasing_jobs = Job()
    job.addChild(split_job)
    split_job.addFollowOn(phasing_jobs)
    
    # Extract contig names from path_list
    work_dir = job.fileStore.getLocalTempDir()
    path_list_file_path = os.path.join(work_dir, os.path.basename(path_list_id))
    job.fileStore.readGlobalFile(path_list_id, path_list_file_path)
    contigs_list = []
    with open(path_list_file_path) as in_path_names:
        for line in in_path_names:
            toks = line.strip().split()
            contig_id = toks[0]
            contigs_list.append(contig_id)
    
    split_jointcalled_vcf_job = split_job.addChildJobFn(run_split_jointcalled_vcf, context, joint_called_vcf_id, joint_called_vcf_index_id, proband_name, maternal_name, paternal_name, contigs_list, filter_parents=False)
    
    phased_vcf_ids = []
    phased_vcf_index_ids = []
    phased_vcf_names = []
    for contig_id in contigs_list:
        phasing_job = phasing_jobs.addChildJobFn(run_whatshap_phasing, context, split_jointcalled_vcf_job.rv(contig_id), contig_id, proband_name, maternal_name, paternal_name,
                                                    proband_bam_id, proband_bam_index_id,
                                                    maternal_bam_id, maternal_bam_index_id,
                                                    paternal_bam_id, paternal_bam_index_id,
                                                    ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                    ped_file_id, genetic_map_id,
                                                    cores=context.config.misc_cores,
                                                    memory=context.config.alignment_mem,
                                                    disk=context.config.alignment_disk)
        phased_vcf_ids.append(phasing_job.rv(0))
        phased_vcf_index_ids.append(phasing_job.rv(1))
        phased_vcf_names.append(phasing_job.rv(2))
    
    construct_job = phasing_jobs.addFollowOnJobFn(run_construct_index_workflow, context, options, '{}.parental.graphs'.format(proband_name), ref_fasta_id,
                                                    phased_vcf_ids, use_haplotypes=True, use_decoys=options.use_decoys, contigs_list=contigs_list)
    return construct_job.rv()

##################################################
########## VG_WDL CONSTRUCT PORT #################
def run_construct_index_workflow(job, context, options, graph_name, ref_fasta_id, contig_vcf_gz_id_list,
                                    use_haplotypes=False, use_decoys=False, 
                                    contigs_list=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y","MT"],
                                    decoy_regex=">GL\|>NC_007605\|>hs37d5"):
    # we make a sub job tree so that all graph construction and indexing is encapsulated in a top-level job
    RealtimeLogger.info("Running run_construct_index_workflow")
    construct_jobs = Job()
    indexing_jobs = Job()
    job.addChild(construct_jobs)
    construct_indexes = {}
    
    construct_chromosome_graph_vg_ids = []
    for contig, vcf_gz_id in zip(contigs_list,contig_vcf_gz_id_list):
        construct_chromosome_graph_vg_ids.append(construct_jobs.addChildJobFn(run_construct_graph_pedigree, context, options,
                                        ref_fasta_id, contig, vcf_gz_id=vcf_gz_id, use_haplotypes=use_haplotypes,
                                        cores=context.config.construct_cores,
                                        memory=context.config.fq_split_mem,
                                        disk=context.config.calling_disk).rv())
    
    construct_decoy_graph_vg_ids = []
    if use_decoys:
        extract_decoys_job = construct_jobs.addChildJobFn(run_extract_decoys, context, options, ref_fasta_id)
        construct_decoy_graph_vg_ids_job = extract_decoys_job.addFollowOnJobFn(run_construct_decoy_contigs_subworkflow, context, options, ref_fasta_id, extract_decoys_job.rv(),
                                                                                 cores=context.config.construct_cores,
                                                                                 memory=context.config.fq_split_mem,
                                                                                 disk=context.config.calling_disk)
        combine_graphs_job = construct_jobs.addFollowOnJobFn(run_combine_graphs, context, options, graph_name, construct_chromosome_graph_vg_ids, decoy_contigs_vg_id_list=construct_decoy_graph_vg_ids_job.rv(),
                                                                cores=context.config.construct_cores,
                                                                memory=context.config.construct_mem,
                                                                disk=context.config.construct_disk)
    else:
        combine_graphs_job = construct_jobs.addFollowOnJobFn(run_combine_graphs, context, options, graph_name, construct_chromosome_graph_vg_ids,
                                                                cores=context.config.construct_cores,
                                                                memory=context.config.construct_mem,
                                                                disk=context.config.construct_disk)
        
    combine_graphs_job.addFollowOn(indexing_jobs)
    xg_index_job = indexing_jobs.addChildJobFn(run_xg_index, context, options, graph_name, combine_graphs_job.rv(0),
                                                    cores=context.config.alignment_cores,
                                                    memory=context.config.construct_mem,
                                                    disk=context.config.construct_disk)
    
    index_output = {}
    construct_indexes = {}
    construct_indexes['vg'] = combine_graphs_job.rv(0)
    construct_indexes['xg'] = xg_index_job.rv()
    if use_haplotypes:
        contig_gbwt_ids = []
        contig_gbwt_ids_job = indexing_jobs.addChildJobFn(run_gbwt_index_subworkflow, context, options, combine_graphs_job.rv(2), contig_vcf_gz_id_list)
        gbwt_merge_job = contig_gbwt_ids_job.addFollowOnJobFn(run_gbwt_merge, context, options, contig_gbwt_ids_job.rv(), graph_name,
                                                                cores=context.config.construct_cores,
                                                                memory=context.config.alignment_mem,
                                                                disk=context.config.construct_disk)
        prune_graph_with_haplotypes_job = gbwt_merge_job.addChildJobFn(run_prune_graph_with_haplotypes, context, options, combine_graphs_job.rv(2), contig_gbwt_ids_job.rv(), combine_graphs_job.rv(1),
                                                                cores=context.config.prune_cores,
                                                                memory=context.config.prune_mem,
                                                                disk=context.config.prune_disk)
        if use_decoys:
            prune_decoy_graph_jobs = gbwt_merge_job.addChildJobFn(run_prune_graph_subworkflow, context, options, combine_graphs_job.rv(4))
            gcsa_index_job = gbwt_merge_job.addFollowOnJobFn(run_gcsa_index, context, options, graph_name, prune_graph_with_haplotypes_job.rv(0), prune_graph_with_haplotypes_job.rv(1), prune_decoy_graph_jobs.rv(),
                                                                cores=context.config.alignment_cores,
                                                                memory=context.config.gcsa_index_mem,
                                                                disk=context.config.gcsa_index_disk)
        else:
            gcsa_index_job = gbwt_merge_job.addFollowOnJobFn(run_gcsa_index, context, options, graph_name, prune_graph_with_haplotypes_job.rv(0), prune_graph_with_haplotypes_job.rv(1),
                                                                cores=context.config.alignment_cores,
                                                                memory=context.config.gcsa_index_mem,
                                                                disk=context.config.gcsa_index_disk)
        construct_indexes['gcsa'] = gcsa_index_job.rv(0)
        construct_indexes['lcp'] = gcsa_index_job.rv(1)
        construct_indexes['gbwt'] = gbwt_merge_job.rv()
    else:
        prune_graph_ids = []
        prune_graph_ids_job = indexing_jobs.addChildJobFn(run_prune_graph_subworkflow, context, options, combine_graphs_job.rv(3))
        gcsa_index_job = prune_graph_ids_job.addFollowOnJobFn(run_gcsa_index, context, options, graph_name, prune_graph_ids_job.rv(), combine_graphs_job.rv(1),
                                                                cores=context.config.alignment_cores,
                                                                memory=context.config.gcsa_index_mem,
                                                                disk=context.config.gcsa_index_disk)
        construct_indexes['gcsa'] = gcsa_index_job.rv(0)
        construct_indexes['lcp'] = gcsa_index_job.rv(1)
    return construct_indexes

#def run_construct_decoy_contigs_subworkflow(job, context, options, ref_fasta_id, decoy_contigs_file_id):
def run_construct_decoy_contigs_subworkflow(job, context, options, ref_fasta_id, decoy_contigs_list):
    child_jobs = Job()
    job.addChild(child_jobs)
    RealtimeLogger.info("Running run_construct_decoy_contigs_subworkflow, decoy_contigs: {}".format(str(decoy_contigs_list)))
    RealtimeLogger.info("Running run_construct_decoy_contigs_subworkflow, ref_fasta_id: {}".format(ref_fasta_id))
    construct_decoy_graph_vg_ids = []
    for decoy_contig in decoy_contigs_list:
        construct_decoy_graph_vg_ids.append(child_jobs.addChildJobFn(run_construct_graph_pedigree, context, options, 
                                            ref_fasta_id, decoy_contig, use_haplotypes=False,
                                            cores=context.config.construct_cores,
                                            memory=context.config.fq_split_mem,
                                            disk=context.config.calling_disk).rv())
    
    return construct_decoy_graph_vg_ids

def run_gbwt_index_subworkflow(job, context, options, vg_ids, contig_vcf_gz_id_list):
    child_jobs = Job()
    job.addChild(child_jobs)
    RealtimeLogger.info("Running run_gbwt_index_subworkflow")
    contig_gbwt_ids = []
    for vg_id, vcf_gz_id in zip(vg_ids, contig_vcf_gz_id_list):
        contig_gbwt_ids.append(child_jobs.addChildJobFn(run_gbwt_index, context, options, vg_id, vcf_gz_id,
                                                                cores=context.config.gbwt_index_cores,
                                                                memory=context.config.gbwt_index_mem,
                                                                disk=context.config.gbwt_index_disk).rv())
    
    return contig_gbwt_ids

def run_prune_graph_subworkflow(job, context, options, contig_vg_ids):
    child_jobs = Job()
    job.addChild(child_jobs)
    RealtimeLogger.info("Running run_prune_graph_subworkflow")
    prune_graph_ids = []
    for contig_vg_id in contig_vg_ids:
        prune_graph_ids.append(child_jobs.addChildJobFn(run_prune_graph, context, options, contig_vg_id,
                                                                cores=context.config.prune_cores,
                                                                memory=context.config.prune_mem,
                                                                disk=context.config.prune_disk).rv())
    
    return prune_graph_ids

def run_extract_decoys(job, context, options, ref_fasta_id, decoy_regex='\'>GL|>NC_007605|>hs37d5\''):
    RealtimeLogger.info("Running run_extract_decoys")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    ref_fasta_path =  os.path.join(work_dir, os.path.basename(ref_fasta_id))
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    cmd = [['cat', os.path.basename(ref_fasta_path)]]
    cmd.append(['grep', '-E', decoy_regex])
    cmd.append(['cut', '-f', '1', '-d', '\' \''])
    cmd.append(['cut', '-f', '2', '-d', '\'>\'', '>>', 'decoy_contig_ids.txt'])
    chain_cmds = [' '.join(p) for p in cmd]
    command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' | '.join(chain_cmds))]
    context.runner.call(job, command, work_dir = work_dir, tool_name='vg')
    decoy_contigs_list = []
    with open(os.path.join(work_dir, 'decoy_contig_ids.txt'), 'r') as output_contigs_list:
        for line in output_contigs_list:
            decoy_contigs_list.append(line.strip())
    
    context.write_output_file(job, os.path.join(work_dir, 'decoy_contig_ids.txt'))
    return decoy_contigs_list

def run_construct_graph_pedigree(job, context, options, ref_fasta_id, contig_name, vcf_gz_id=None, use_haplotypes=False, vg_construct_options="--node-max 32 --handle-sv"):
    RealtimeLogger.info("Running run_construct_graph_pedigree")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    ref_fasta_path = os.path.join(work_dir, os.path.basename(ref_fasta_id))
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    command = ['vg', 'construct', '--threads', job.cores, '-R', contig_name, '-C', '-r', os.path.basename(ref_fasta_path)]
    if vcf_gz_id is not None:
        vcf_gz_path = os.path.join(work_dir, os.path.basename(vcf_gz_id))
        job.fileStore.readGlobalFile(vcf_gz_id, vcf_gz_path)
        context.runner.call(job, ['tabix', os.path.basename(vcf_gz_path)], work_dir = work_dir, tool_name='vg')
        command += ['-v', os.path.basename(vcf_gz_path), '--region-is-chrom']
        
    vg_construct_options_list = vg_construct_options.split()
    command += vg_construct_options.split()
    if use_haplotypes:
        command += ['-a']
    with open(os.path.join(work_dir, '{}.vg'.format(contig_name)), 'wb') as output_contig_vg_file:
        context.runner.call(job, command, work_dir = work_dir, tool_name='vg', outfile=output_contig_vg_file)
    context.runner.call(job, ['rm', '-f', os.path.basename(ref_fasta_path)], work_dir = work_dir, tool_name='vg')
    
    # Write output to the outstore
    contig_vg_id = context.write_intermediate_file(job, os.path.join(work_dir, '{}.vg'.format(contig_name)))
    return contig_vg_id

def run_combine_graphs(job, context, options, graph_name, contigs_vg_id_list, decoy_contigs_vg_id_list=[]):
    RealtimeLogger.info("Running run_combine_graphs")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    contigs_vg_id_list+=decoy_contigs_vg_id_list
    decoy_contigs_uid_vg = []
    contigs_uid_vg = []
    all_contigs_uid_vg = []
    all_contigs_uid_vg_ids = []
    for contig_vg_id in contigs_vg_id_list:
        contig_vg_path = os.path.join(work_dir, os.path.basename(contig_vg_id))
        job.fileStore.readGlobalFile(contig_vg_id, contig_vg_path)
        contig_name = os.path.splitext(os.path.basename(contig_vg_path))[0]
        if ("GL" in contig_name) or ("NC_007605" in contig_name) or ("hs37d5" in contig_name):
            decoy_contigs_uid_vg.append(contig_vg_id)
        else:
            contigs_uid_vg.append(contig_vg_id)
        all_contigs_uid_vg.append(os.path.basename(contig_vg_path))
        all_contigs_uid_vg_ids.append(contig_vg_id)
    
    command = ['vg', 'ids', '-j', '-m', 'empty.id_map']
    command += all_contigs_uid_vg
    context.runner.call(job, command, work_dir = work_dir, tool_name='vg')
    command = ['vg', 'combine']
    command += all_contigs_uid_vg
    with open(os.path.join(work_dir, '{}.vg'.format(graph_name)), 'wb') as output_merged_vg_file:
        context.runner.call(job, command, work_dir = work_dir, tool_name='vg', outfile=output_merged_vg_file)
    
    combined_vg_id = context.write_output_file(job, os.path.join(work_dir, '{}.vg'.format(graph_name)))
    empty_id_map_id = context.write_output_file(job, os.path.join(work_dir, 'empty.id_map'))
    
    return (combined_vg_id, empty_id_map_id, contigs_uid_vg, all_contigs_uid_vg_ids, decoy_contigs_uid_vg)

def run_gbwt_index(job, context, options, vg_id, vcf_gz_id):
    RealtimeLogger.info("Running run_gbwt_index")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vg_file_path = os.path.join(work_dir, os.path.basename(vg_id))
    job.fileStore.readGlobalFile(vg_id, vg_file_path)
    
    vcf_gz_path = os.path.join(work_dir, os.path.basename(vcf_gz_id))
    job.fileStore.readGlobalFile(vcf_gz_id, vcf_gz_path)
    
    contig_name = os.path.splitext(os.path.basename(vg_file_path))[0]
    context.runner.call(job, ['tabix', os.path.basename(vcf_gz_path)], work_dir = work_dir, tool_name='vg')
    context.runner.call(job, ['vg', 'index', '--threads', job.cores, '--force-phasing', '--discard-overlaps', '-G', '{}.gbwt'.format(contig_name), '-v', os.path.basename(vcf_gz_path), os.path.basename(vg_file_path)], work_dir = work_dir, tool_name='vg')
    
    gbwt_id = context.write_intermediate_file(job, os.path.join(work_dir, '{}.gbwt'.format(contig_name)))
    
    return gbwt_id

def run_gbwt_merge(job, context, options, gbwt_id_list, graph_name):
    RealtimeLogger.info("Running run_gbwt_merge")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    gbwt_list = []
    for gbwt_id in gbwt_id_list:
        gbwt_file_path = os.path.join(work_dir, os.path.basename(gbwt_id))
        job.fileStore.readGlobalFile(gbwt_id, gbwt_file_path)
        gbwt_list.append(os.path.basename(gbwt_file_path))
    
    RealtimeLogger.info("Merging gbwt list: {}".format(str(gbwt_list)))
    outputfile = os.path.join(work_dir, '{}.gbwt'.format(graph_name))
    if len(gbwt_list) > 1:
        command = ['vg', 'gbwt', '-m', '-f', '-o', '{}.gbwt'.format(graph_name)]
        command += gbwt_list
        context.runner.call(job, command, work_dir = work_dir, tool_name='vg')
    else:
        outputfile = os.path.join(work_dir, gbwt_list[0])
    
    merged_gbwt_id = context.write_output_file(job, outputfile)
    
    return merged_gbwt_id

def run_xg_index(job, context, options, graph_name, vg_id, xg_options=None):
    RealtimeLogger.info("Running run_xg_index")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vg_file_path = os.path.join(work_dir, os.path.basename(vg_id))
    job.fileStore.readGlobalFile(vg_id, vg_file_path)
    
    command = ['vg', 'index', '--threads', job.cores, '-x', '{}.xg'.format(graph_name)]
    if xg_options:
        command += [xg_options]
    command += [os.path.basename(vg_file_path)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='vg')
    
    xg_id = context.write_output_file(job, os.path.join(work_dir, '{}.xg'.format(graph_name)))
    
    return xg_id

def run_prune_graph(job, context, options, contig_vg_id, prune_options=None):
    RealtimeLogger.info("Running run_prune_graph")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vg_file_path = os.path.join(work_dir, os.path.basename(contig_vg_id))
    job.fileStore.readGlobalFile(contig_vg_id, vg_file_path)
    contig_name = os.path.splitext(os.path.basename(vg_file_path))[0]
    
    command = ['vg', 'prune', '--threads', job.cores, '-r', os.path.basename(vg_file_path)]
    if prune_options:
        command += [prune_options]
    with open(os.path.join(work_dir, '{}.pruned.vg'.format(contig_name)), 'wb') as output_pruned_vg_file:
        context.runner.call(job, command, work_dir = work_dir, tool_name='vg', outfile=output_pruned_vg_file)
    pruned_vg_id = context.write_intermediate_file(job, os.path.join(work_dir, '{}.pruned.vg'.format(contig_name)))
    
    return pruned_vg_id

def run_prune_graph_with_haplotypes(job, context, options, contig_vg_ids_list, contig_gbwt_ids_list, empty_id_map_id, prune_options=None):
    RealtimeLogger.info("Running run_prune_graph_with_haplotypes")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    empty_id_map_path = os.path.join(work_dir, os.path.basename(empty_id_map_id))
    job.fileStore.readGlobalFile(empty_id_map_id, empty_id_map_path)
    
    pruned_vg_id_list = []
    for vg_id, gbwt_id in zip(contig_vg_ids_list, contig_gbwt_ids_list):
        vg_file_path = os.path.join(work_dir, os.path.basename(vg_id))
        job.fileStore.readGlobalFile(vg_id, vg_file_path)
        contig_name = os.path.splitext(os.path.basename(vg_file_path))[0]
        gbwt_file_path = os.path.join(work_dir, os.path.basename(gbwt_id))
        job.fileStore.readGlobalFile(gbwt_id,  gbwt_file_path)
        command = ['vg', 'prune', '--threads', job.cores, '-u', '-g', os.path.basename(gbwt_file_path), '-a', '-m', os.path.basename(empty_id_map_path), os.path.basename(vg_file_path)]
        if prune_options:
            command += [prune_options]
        with open(os.path.join(work_dir, '{}.pruned.vg'.format(contig_name)), 'wb') as output_pruned_vg_file:
            context.runner.call(job, command, work_dir = work_dir, tool_name='vg', outfile=output_pruned_vg_file)
        pruned_vg_id_list.append(context.write_intermediate_file(job, os.path.join(work_dir, '{}.pruned.vg'.format(contig_name))))
    
    pruned_id_map_id = context.write_output_file(job, os.path.join(work_dir, os.path.basename(empty_id_map_path)))
    
    return (pruned_vg_id_list, pruned_id_map_id)

def run_gcsa_index(job, context, options, graph_name, pruned_vg_ids_list, id_map_id, prune_decoy_vg_ids_list=[], gcsa_options=None):
    RealtimeLogger.info("Running run_gcsa_index")
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    id_map_path = os.path.join(work_dir, os.path.basename(id_map_id))
    job.fileStore.readGlobalFile(id_map_id, id_map_path)
    
    pruned_vg_path_list = []
    pruned_vg_ids_list+=prune_decoy_vg_ids_list
    for pruned_vg_id in pruned_vg_ids_list:
        pruned_vg_path = os.path.join(work_dir, os.path.basename(pruned_vg_id))
        job.fileStore.readGlobalFile(pruned_vg_id, pruned_vg_path)
        pruned_vg_path_list.append(os.path.basename(pruned_vg_path))
    
    command = ['vg', 'index', '--threads', job.cores, '-p', '-g', '{}.gcsa'.format(graph_name), '-f', os.path.basename(id_map_path)]
    if gcsa_options:
        command += [gcsa_options]
    command += pruned_vg_path_list
    context.runner.call(job, command, work_dir = work_dir, tool_name='vg')
    
    gcsa_file_id = context.write_output_file(job, os.path.join(work_dir, '{}.gcsa'.format(graph_name)))
    gcsa_lcp_file_id = context.write_output_file(job, os.path.join(work_dir, '{}.gcsa.lcp'.format(graph_name)))
    
    return (gcsa_file_id, gcsa_lcp_file_id)

##################################################
##################################################
##################################################

def run_process_parental_graph_index(job, context, options, parental_indexes, old_indexes):
    if 'id_ranges' not in parental_indexes.keys():
        parental_indexes['id_ranges'] = old_indexes['id_ranges']
    # Remove the gbwt index for now as it's still causing some mapping performance issues
    parental_indexes.pop('gbwt', None)
    return parental_indexes

def run_snpEff_annotation(job, context, cohort_name, joint_called_vcf_id, snpeff_database_file_id):
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    joint_called_vcf_path =  os.path.join(work_dir, os.path.basename(joint_called_vcf_id))
    job.fileStore.readGlobalFile(joint_called_vcf_id, joint_called_vcf_path)
    
    snpeff_database_file_path =  os.path.join(work_dir, os.path.basename(snpeff_database_file_id))
    job.fileStore.readGlobalFile(snpeff_database_file_id, snpeff_database_file_path)
    
    command = ['bcftools', 'norm', '-m-both', '--threads', job.cores,
                    '-o', '{}.unrolled.vcf'.format(cohort_name), os.path.basename(joint_called_vcf_path)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='bcftools')
    command = ['unzip', os.path.basename(snpeff_database_file_path)]
    context.runner.call(job, command, work_dir = work_dir)
    command = ['snpEff', '-Xmx{}g'.format(int(float(job.memory)/1000000000)), '-i', 'VCF', '-o', 'VCF', '-noLof', '-noHgvs', '-formatEff', '-classic',
                    '-dataDir', '$PWD/data',
                    'GRCh37.75', '{}.unrolled.vcf'.format(cohort_name)]
    with open(os.path.join(work_dir, '{}.snpeff.unrolled.vcf'.format(cohort_name)), 'wb') as output_snpeff_vcf:
        context.runner.call(job, command, work_dir = work_dir, tool_name='snpEff', outfile=output_snpeff_vcf)
    context.runner.call(job, ['bgzip', '{}.snpeff.unrolled.vcf'.format(cohort_name)], work_dir = work_dir, tool_name='vg')
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', '{}.snpeff.unrolled.vcf.gz'.format(cohort_name)], work_dir=work_dir)
    
    # checkpoint to out store
    out_file = os.path.join(work_dir, '{}.snpeff.unrolled.vcf.gz'.format(cohort_name))
    snpeff_annotated_vcf_file_id = context.write_output_file(job, out_file)
    snpeff_annotated_vcf_index_file_id = context.write_output_file(job, out_file + '.tbi')
    
    return (snpeff_annotated_vcf_file_id, snpeff_annotated_vcf_index_file_id)

def run_indel_realignment(job, context, sample_name, sample_bam_id, ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id):
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    sample_bam_path = os.path.join(work_dir, os.path.basename(sample_bam_id))
    job.fileStore.readGlobalFile(sample_bam_id, sample_bam_path)
    
    ref_fasta_name = os.path.basename(ref_fasta_id)
    ref_fasta_path = os.path.join(work_dir, ref_fasta_name)
    job.fileStore.readGlobalFile(ref_fasta_id, ref_fasta_path)
    
    ref_fasta_index_path = os.path.join(work_dir, '{}.fai'.format(ref_fasta_name))
    job.fileStore.readGlobalFile(ref_fasta_index_id, ref_fasta_index_path)
    
    ref_fasta_dict_path = os.path.join(work_dir, '{}.dict'.format(os.path.splitext(ref_fasta_name)[0]))
    job.fileStore.readGlobalFile(ref_fasta_dict_id, ref_fasta_dict_path)
    
    
    command = ['samtools', 'index', os.path.basename(sample_bam_path)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='samtools')
    command = ['java', '-jar', '/usr/GenomeAnalysisTK.jar', '-T', 'RealignerTargetCreator',
               '--remove_program_records', '-drf', 'DuplicateRead', '--disable_bam_indexing',
               '-nt', str(job.cores), '-R', os.path.basename(ref_fasta_path), '-I', os.path.basename(sample_bam_path),
               '--out', 'forIndelRealigner.intervals']
    context.runner.call(job, command, work_dir = work_dir, tool_name='gatk3')
    command = ['java', '-jar', '/usr/GenomeAnalysisTK.jar', '-T', 'IndelRealigner',
               '--remove_program_records', '-drf', 'DuplicateRead', '--disable_bam_indexing',
               '-R', os.path.basename(ref_fasta_path), '-I', os.path.basename(sample_bam_path),
               '--targetIntervals', 'forIndelRealigner.intervals',
               '--out', '{}.indel_realigned.bam'.format(sample_name)]
    context.runner.call(job, command, work_dir = work_dir, tool_name='gatk3')
    cmd_list = []
    cmd_list.append(['samtools', 'sort', '--threads', str(job.cores), '-O', 'BAM',
               '{}.indel_realigned.bam'.format(sample_name)])
    cmd_list.append(['samtools', 'calmd', '-b', '-', os.path.basename(ref_fasta_path)])
    with open(os.path.join(work_dir, '{}_positionsorted.mdtag.indel_realigned.bam'.format(os.path.basename(sample_bam_id))), 'wb') as output_indel_realigned_bam:
        context.runner.call(job, cmd_list, work_dir = work_dir, tool_name='samtools', outfile=output_indel_realigned_bam)
    
    # Delete input files
    job.fileStore.deleteGlobalFile(sample_bam_id)
    
    return context.write_intermediate_file(job, os.path.join(work_dir, '{}_positionsorted.mdtag.indel_realigned.bam'.format(os.path.basename(sample_bam_id))))

def run_cohort_indel_realign_pipeline(job, context, options, proband_name, maternal_name, paternal_name,
                                        proband_chr_bam_ids, maternal_chr_bam_ids, paternal_chr_bam_ids,
                                        path_list, ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                        siblings_names=None, sibling_mapping_chr_bam_ids_dict=None):
    
    # we make a child jobs so that all indel realignment is encapsulated in a top-level job
    RealtimeLogger.info("Starting cohort indel realignment workflow")
    proband_indel_realign_jobs = Job()
    job.addChild(proband_indel_realign_jobs)
    maternal_indel_realign_jobs = Job()
    job.addChild(maternal_indel_realign_jobs)
    paternal_indel_realign_jobs = Job()
    job.addChild(paternal_indel_realign_jobs)
    
    sibling_merge_bam_ids_list = None
    sibling_merge_bam_index_ids_list = None
    if siblings_names is not None:
        sibling_root_job_dict =  {}
        sibling_indel_realign_job_dict = {}
        sibling_merge_bam_ids_list = []
        for sibling_name in siblings_names:
            # Dynamically allocate sibling indel realignment jobs to overall workflow structure
            sibling_root_job_dict[sibling_name] = Job()
            job.addChild(sibling_root_job_dict[sibling_name])
            sibling_chr_bam_indel_realign_output = []
            for sibling_chr_bam_id in sibling_mapping_chr_bam_ids_dict[sibling_name]:
                sibling_indel_realign_job_dict[sibling_name] = sibling_root_job_dict[sibling_name].addChildJobFn(run_indel_realignment, context,
                                                                        sibling_name, sibling_chr_bam_id,
                                                                        ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                                        cores=context.config.alignment_cores,
                                                                        memory=context.config.alignment_mem,
                                                                        disk=context.config.alignment_disk)
                sibling_chr_bam_indel_realign_output.append(sibling_indel_realign_job_dict[sibling_name].rv())
            sibling_merge_bam_ids_list.append(sibling_root_job_dict[sibling_name].addFollowOnJobFn(run_merge_bams_ped_workflow, context, sibling_name, sibling_chr_bam_indel_realign_output, True, True,
                                                                                                                cores=context.config.alignment_cores, memory=context.config.alignment_mem, disk=context.config.alignment_disk).rv())
    
    proband_chr_bam_indel_realign_output = []
    for proband_chr_bam_id in proband_chr_bam_ids:
        proband_chr_bam_indel_realign_job = proband_indel_realign_jobs.addChildJobFn(run_indel_realignment, context,
                                                                proband_name, proband_chr_bam_id,
                                                                ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                                cores=context.config.alignment_cores,
                                                                memory=context.config.alignment_mem,
                                                                disk=context.config.alignment_disk)
        proband_chr_bam_indel_realign_output.append(proband_chr_bam_indel_realign_job.rv())
    proband_merged_bam_job = proband_indel_realign_jobs.addFollowOnJobFn(run_merge_bams_ped_workflow, context, proband_name, proband_chr_bam_indel_realign_output, True, True,
                                                                                        cores=context.config.alignment_cores, memory=context.config.alignment_mem, disk=context.config.alignment_disk)
    
    maternal_chr_bam_indel_realign_output = []
    for maternal_chr_bam_id in maternal_chr_bam_ids:
        maternal_chr_bam_indel_realign_job = maternal_indel_realign_jobs.addChildJobFn(run_indel_realignment, context,
                                                                maternal_name, maternal_chr_bam_id,
                                                                ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                                cores=context.config.alignment_cores,
                                                                memory=context.config.alignment_mem,
                                                                disk=context.config.alignment_disk)
        maternal_chr_bam_indel_realign_output.append(maternal_chr_bam_indel_realign_job.rv())
    maternal_merged_bam_job = maternal_indel_realign_jobs.addFollowOnJobFn(run_merge_bams_ped_workflow, context, maternal_name, maternal_chr_bam_indel_realign_output, True, True,
                                                                                        cores=context.config.alignment_cores, memory=context.config.alignment_mem, disk=context.config.alignment_disk)
    
    paternal_chr_bam_indel_realign_output = []
    for paternal_chr_bam_id in paternal_chr_bam_ids:
        paternal_chr_bam_indel_realign_job = paternal_indel_realign_jobs.addChildJobFn(run_indel_realignment, context,
                                                                paternal_name, paternal_chr_bam_id,
                                                                ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id,
                                                                cores=context.config.alignment_cores,
                                                                memory=context.config.alignment_mem,
                                                                disk=context.config.alignment_disk)
        paternal_chr_bam_indel_realign_output.append(paternal_chr_bam_indel_realign_job.rv())
    paternal_merged_bam_job = paternal_indel_realign_jobs.addFollowOnJobFn(run_merge_bams_ped_workflow, context, paternal_name, paternal_chr_bam_indel_realign_output, True, True,
                                                                                        cores=context.config.alignment_cores, memory=context.config.alignment_mem, disk=context.config.alignment_disk)
    
    return (proband_merged_bam_job.rv(0), proband_merged_bam_job.rv(1), maternal_merged_bam_job.rv(0), maternal_merged_bam_job.rv(1), paternal_merged_bam_job.rv(0), paternal_merged_bam_job.rv(1), sibling_merge_bam_ids_list)

def run_collect_sibling_indel_realigned_bams(job, context, options, sibling_merge_bam_bamindex_id_pairs):
    RealtimeLogger.debug("Collecting sibling bams: {}".format(sibling_merge_bam_bamindex_id_pairs))
    sibling_bam_ids_list = []
    sibling_bam_index_ids_list = []
    for bam_bam_index_pair in sibling_merge_bam_bamindex_id_pairs:
        sibling_bam_ids_list.append(bam_bam_index_pair[0])
        sibling_bam_index_ids_list.append(bam_bam_index_pair[1])

    return (sibling_bam_ids_list, sibling_bam_index_ids_list)
    

def run_pedigree(job, context, options, fastq_proband, gam_input_reads_proband, bam_input_reads_proband,
                fastq_maternal, gam_input_reads_maternal, bam_input_reads_maternal,
                fastq_paternal, gam_input_reads_paternal, bam_input_reads_paternal,
                fastq_siblings, gam_input_reads_siblings, bam_input_reads_siblings,
                proband_name, maternal_name, paternal_name, siblings_names,
                indel_realign_bams, snpeff_annotation, interleaved, mapper, indexes, ref_fasta_ids, misc_file_ids,
                reads_file_ids_proband=None, reads_file_ids_maternal=None, reads_file_ids_paternal=None, reads_file_ids_siblings=None,
                reads_chunk_ids_proband=None, reads_chunk_ids_maternal=None, reads_chunk_ids_paternal=None, reads_chunk_ids_siblings=None,
                bam_output=False, surject=False, 
                gbwt_penalty=None, validate=False):
    """
    Split the fastq, then align each chunk.
    
    Exactly one of fastq, gam_input_reads, or bam_input_reads should be
    non-falsey, to indicate what kind of data the file IDs in reads_file_ids or
    reads_chunk_ids correspond to.
    
    Exactly one of reads_file_ids or read_chunks_ids should be specified.
    reads_file_ids holds a list of file IDs of non-chunked input read files,
    which will be chunked if necessary. reads_chunk_ids holds lists of chunk
    IDs for each read file, as produced by run_split_reads_if_needed.
    
    indexes is a dict from index type ('xg', 'gcsa', 'lcp', 'id_ranges',
    'gbwt', 'minimizer', 'distance', 'snarls') to index file ID. Some indexes
    are extra and specifying them will change mapping behavior. Some indexes
    are required for certain values of mapper.
    
    mapper can be 'map', 'mpmap', or 'gaffe'. For 'map' and 'mpmap', the 'gcsa'
    and 'lcp' indexes are required. For 'gaffe', the 'gbwt', 'minimizer' and
    'distance' indexes are required. All the mappers require the 'xg' index.
    
    If bam_output is set, produce BAMs. If surject is set, surject reads down
    to paths. 
    
    If the 'gbwt' index is present and gbwt_penalty is specified, the default
    recombination penalty will be overridden.
    
    returns output gams, one per chromosome, the total mapping time (excluding
    toil-vg overhead such as transferring and splitting files), and output
    BAMs, one per chromosome, if computed.
    """
    
    # Make sure we have exactly one type of input
    assert (bool(fastq_proband) + bool(gam_input_reads_proband) + bool(bam_input_reads_proband) == 1)
    assert (bool(fastq_maternal) + bool(gam_input_reads_maternal) + bool(bam_input_reads_maternal) == 1)
    assert (bool(fastq_paternal) + bool(gam_input_reads_paternal) + bool(bam_input_reads_paternal) == 1)
    assert (bool(fastq_siblings) + bool(gam_input_reads_siblings) + bool(bam_input_reads_siblings) <= 1)
    
    # Make sure we have exactly one kind of file IDs
    assert(bool(reads_file_ids_proband) + bool(reads_chunk_ids_proband) == 1)
    assert(bool(reads_file_ids_maternal) + bool(reads_chunk_ids_maternal) == 1)
    assert(bool(reads_file_ids_paternal) + bool(reads_chunk_ids_paternal) == 1)
    assert(bool(reads_file_ids_siblings) + bool(reads_chunk_ids_siblings) <= 1)
    
    # define the job workflow structure
    stage1_jobs = Job()
    stage2_jobs = Job()
    stage3_jobs = Job()
    stage4_jobs = Job()
    
    job.addChild(stage1_jobs)
    stage1_jobs.addFollowOn(stage2_jobs)
    stage2_jobs.addFollowOn(stage3_jobs)
    stage3_jobs.addFollowOn(stage4_jobs)
    
    
    # Define the probands 1st alignment and variant calling jobs
    proband_first_mapping_job = stage1_jobs.addChildJobFn(run_mapping, context, fastq_proband,
                                     gam_input_reads_proband, bam_input_reads_proband,
                                     proband_name,
                                     interleaved, mapper, indexes,
                                     reads_file_ids_proband,
                                     bam_output=bam_output, surject=surject,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    proband_calling_job = proband_first_mapping_job.addFollowOnJobFn(run_pipeline_call_gvcfs, context, options,
                                proband_name, proband_first_mapping_job.rv(2), 
                                ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                run_dragen=options.run_dragen,
                                dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                cores=context.config.misc_cores,
                                memory=context.config.misc_mem,
                                disk=context.config.misc_disk)
    
    # Define the maternal alignment and variant calling jobs
    maternal_mapping_job = stage1_jobs.addChildJobFn(run_mapping, context, fastq_maternal,
                                     gam_input_reads_maternal, bam_input_reads_maternal,
                                     maternal_name,
                                     interleaved, mapper, indexes,
                                     reads_file_ids_maternal,
                                     bam_output=options.bam_output, surject=options.surject,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    maternal_calling_job = maternal_mapping_job.addFollowOnJobFn(run_pipeline_call_gvcfs, context, options,
                                maternal_name, maternal_mapping_job.rv(2),
                                ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                run_dragen=options.run_dragen,
                                dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                cores=context.config.misc_cores,
                                memory=context.config.misc_mem,
                                disk=context.config.misc_disk)
    # Dragen calling depends on previous call execution
    #if options.run_dragen:
    #    proband_calling_job.addFollowOn(maternal_calling_job)
    
    # Define the paternal alignment and variant calling jobs
    paternal_mapping_job = stage1_jobs.addChildJobFn(run_mapping, context, fastq_paternal,
                                     gam_input_reads_paternal, bam_input_reads_paternal,
                                     paternal_name,
                                     interleaved, mapper, indexes,
                                     reads_file_ids_paternal,
                                     bam_output=options.bam_output, surject=options.surject,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    paternal_calling_job = paternal_mapping_job.addFollowOnJobFn(run_pipeline_call_gvcfs, context, options,
                                paternal_name, paternal_mapping_job.rv(2),
                                ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                run_dragen=options.run_dragen,
                                dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                cores=context.config.misc_cores,
                                memory=context.config.misc_mem,
                                disk=context.config.misc_disk)
    # Dragen calling depends on previous call execution
    #if options.run_dragen:
    #    maternal_calling_job.addFollowOn(paternal_calling_job)
    
    joint_calling_job = stage2_jobs.addChildJobFn(run_joint_genotyper, context, proband_name,
                                proband_calling_job.rv(2), proband_calling_job.rv(3),
                                maternal_calling_job.rv(2), maternal_calling_job.rv(3),
                                paternal_calling_job.rv(2), paternal_calling_job.rv(3),
                                None, None,
                                ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                snpeff_annotation=snpeff_annotation,
                                run_dragen=options.run_dragen,
                                dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                cores=context.config.calling_cores,
                                memory=context.config.calling_mem,
                                disk=context.config.calling_disk)
    
    gen_map_id = None
    if options.genetic_map is not None:
        gen_map_id = misc_file_ids['genetic_map']
    graph_construction_job = joint_calling_job.addFollowOnJobFn(run_pipeline_construct_parental_graphs, context, options,
                                joint_calling_job.rv(0), joint_calling_job.rv(1),
                                proband_name, maternal_name, paternal_name,
                                proband_calling_job.rv(0), proband_calling_job.rv(1),
                                maternal_calling_job.rv(0), maternal_calling_job.rv(1),
                                paternal_calling_job.rv(0), paternal_calling_job.rv(1),
                                ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                misc_file_ids['path_list'], misc_file_ids['ped_file'], gen_map_id)
    
    # Make a parental graph index collection
    #process_parental_graph_indexes_job = graph_construction_job.addFollowOnJobFn(run_process_parental_graph_index, context, options, graph_construction_job.rv(0,2), indexes)
    process_parental_graph_indexes_job = graph_construction_job.addFollowOnJobFn(run_process_parental_graph_index, context, options, graph_construction_job.rv(), indexes)
    
    # Run stage3 asynchronously for each sample
    # Define the probands 2nd alignment and variant calling jobs
    proband_second_mapping_job = stage3_jobs.addChildJobFn(run_mapping, context, fastq_proband,
                                     gam_input_reads_proband, bam_input_reads_proband,
                                     proband_name,
                                     interleaved, mapper, process_parental_graph_indexes_job.rv(),
                                     reads_file_ids_proband,
                                     bam_output=options.bam_output, surject=options.surject,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    proband_parental_calling_job = proband_second_mapping_job.addFollowOnJobFn(run_pipeline_call_gvcfs, context, options,
                                    proband_name, proband_second_mapping_job.rv(2), 
                                    ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                    run_dragen=options.run_dragen,
                                    dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                    cores=context.config.misc_cores,
                                    memory=context.config.misc_mem,
                                    disk=context.config.misc_disk)
    
    # Define the siblings alignment and variant calling jobs
    sibling_mapping_chr_bam_ids_dict = None
    sibling_merged_bam_ids = None
    sibling_merged_bam_index_ids = None
    sibling_call_gvcf_ids = None
    sibling_call_gvcf_index_ids = None
    if siblings_names is not None: 
        sibling_mapping_job_dict =  {}
        sibling_calling_job_dict =  {}
        sibling_mapping_chr_bam_ids_dict = {}
        sibling_merged_bam_ids = []
        sibling_merged_bam_index_ids = []
        sibling_call_gvcf_ids = []
        sibling_call_gvcf_index_ids = []
        for sibling_number,sibling_name in enumerate(siblings_names):
            
            # Dynamically allocate sibling map allignment jobs to overall workflow structure
            fastq_siblings_collection = None
            gam_input_reads_siblings_collection = None
            bam_input_reads_siblings_collection = None
            reads_file_ids_siblings_list = []
            if fastq_siblings:
                fastq_siblings_collection = fastq_siblings[sibling_number]
                reads_file_ids_siblings_list.append(reads_file_ids_siblings[sibling_number*2:(sibling_number*2)+2])
            elif gam_input_reads_siblings or bam_input_reads_siblings:
                if gam_input_reads_siblings: gam_input_reads_siblings_collection = gam_input_reads_siblings[sibling_number]
                if bam_input_reads_siblings: bam_input_reads_siblings_collection = bam_input_reads_siblings[sibling_number]
                reads_file_ids_siblings_list.append(reads_file_ids_siblings[sibling_number])
            sibling_mapping_job_dict[sibling_name] = stage3_jobs.addChildJobFn(run_mapping, context, fastq_siblings_collection,
                                             gam_input_reads_siblings_collection, bam_input_reads_siblings_collection,
                                             siblings_names[sibling_number],
                                             options.interleaved, options.mapper, process_parental_graph_indexes_job.rv(),
                                             reads_file_ids_siblings_list,
                                             bam_output=options.bam_output, surject=options.surject,
                                             cores=context.config.misc_cores,
                                             memory=context.config.misc_mem,
                                             disk=context.config.misc_disk)
            sibling_calling_job_dict[sibling_name] = sibling_mapping_job_dict[sibling_name].addFollowOnJobFn(run_pipeline_call_gvcfs, context, options,
                                            sibling_name, sibling_mapping_job_dict[sibling_name].rv(2),
                                            ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                            run_dragen=options.run_dragen,
                                            dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                            cores=context.config.misc_cores,
                                            memory=context.config.misc_mem,
                                            disk=context.config.misc_disk)
            # Dragen calling depends on previous call execution
            #if options.run_dragen:
            #    if sibling_number == 0:
            #        proband_parental_calling_job.addFollowOn(sibling_calling_job_dict[sibling_name])
            #    else:
            #        previous_sibling_name = siblings_names[sibling_number-1]
            #        sibling_calling_job_dict[previous_sibling_name].addFollowOn(sibling_calling_job_dict[sibling_name])
            
            # Extract outputs
            sibling_mapping_chr_bam_ids_dict[sibling_name] = sibling_calling_job_dict[sibling_name].rv(4)
            sibling_merged_bam_ids.append(sibling_calling_job_dict[sibling_name].rv(0))
            sibling_merged_bam_index_ids.append(sibling_calling_job_dict[sibling_name].rv(1))
            sibling_call_gvcf_ids.append(sibling_calling_job_dict[sibling_name].rv(2))
            sibling_call_gvcf_index_ids.append(sibling_calling_job_dict[sibling_name].rv(3))
    
    pedigree_joint_call_job = stage4_jobs.addChildJobFn(run_joint_genotyper, context, proband_name,
                                        proband_parental_calling_job.rv(2), proband_parental_calling_job.rv(3),
                                        maternal_calling_job.rv(2), maternal_calling_job.rv(3),
                                        paternal_calling_job.rv(2), paternal_calling_job.rv(3),
                                        sibling_call_gvcf_ids, sibling_call_gvcf_index_ids,
                                        ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                        run_dragen=options.run_dragen,
                                        dragen_ref_index_name=options.dragen_ref_index_name, udp_data_dir=options.udp_data_dir, helix_username=options.helix_username,
                                        snpeff_annotation=snpeff_annotation,
                                        cores=context.config.calling_cores,
                                        memory=context.config.calling_mem,
                                        disk=context.config.calling_disk)
    
    if indel_realign_bams:
        indel_realign_job = stage4_jobs.addChildJobFn(run_cohort_indel_realign_pipeline, context, options,
                                                                proband_name, maternal_name, paternal_name,
                                                                proband_parental_calling_job.rv(4), maternal_calling_job.rv(4), paternal_calling_job.rv(4),
                                                                misc_file_ids['path_list'], ref_fasta_ids[0], ref_fasta_ids[1], ref_fasta_ids[2],
                                                                siblings_names=siblings_names, sibling_mapping_chr_bam_ids_dict=sibling_mapping_chr_bam_ids_dict)
    if snpeff_annotation:
        snpeff_job = pedigree_joint_call_job.addFollowOnJobFn(run_snpEff_annotation, context,
                                                        proband_name, pedigree_joint_call_job.rv(0),
                                                        misc_file_ids['snpeff_data'],
                                                        cores=context.config.calling_cores,
                                                        memory=context.config.calling_mem,
                                                        disk=context.config.calling_disk)
    # Test outputs
    trio_joint_called_vcf = joint_calling_job.rv(0)
    trio_joint_called_vcf_index = joint_calling_job.rv(1)
    # Collect final outputs
    final_proband_bam = proband_parental_calling_job.rv(0)
    final_proband_bam_index = proband_parental_calling_job.rv(1)
    final_proband_gvcf = proband_parental_calling_job.rv(2)
    final_proband_gvcf_index = proband_parental_calling_job.rv(3)
    
    final_maternal_bam = maternal_calling_job.rv(0)
    final_maternal_bam_index = maternal_calling_job.rv(1)
    final_maternal_gvcf = maternal_calling_job.rv(2)
    final_maternal_gvcf_index = maternal_calling_job.rv(3)
    
    final_paternal_bam = paternal_calling_job.rv(0)
    final_paternal_bam_index = paternal_calling_job.rv(1)
    final_paternal_gvcf = paternal_calling_job.rv(2)
    final_paternal_gvcf_index = paternal_calling_job.rv(3)
    
    final_pedigree_joint_called_vcf = pedigree_joint_call_job.rv(0)
    final_pedigree_joint_called_vcf_index = pedigree_joint_call_job.rv(1)
    if snpeff_annotation:
        final_pedigree_joint_called_vcf = snpeff_job.rv(0)
        final_pedigree_joint_called_vcf_index = snpeff_job.rv(1)
    
    final_sibling_bam_list = None
    final_sibling_bam_index_list = None
    final_sibling_gvcf_list = None
    final_sibling_gvcf_index_list = None
    if siblings_names is not None:
        final_sibling_bam_list = sibling_merged_bam_ids
        final_sibling_bam_index_list = sibling_merged_bam_index_ids
        final_sibling_gvcf_list = sibling_call_gvcf_ids
        final_sibling_gvcf_index_list = sibling_call_gvcf_index_ids
    
    if indel_realign_bams:
        final_proband_bam = indel_realign_job.rv(0)
        final_proband_bam_index = indel_realign_job.rv(1)
        final_maternal_bam = indel_realign_job.rv(2)
        final_maternal_bam_index = indel_realign_job.rv(3)
        final_paternal_bam = indel_realign_job.rv(4)
        final_paternal_bam_index = indel_realign_job.rv(5)
        final_sibling_bam_list_job = indel_realign_job.addFollowOnJobFn(run_collect_sibling_indel_realigned_bams, context, options,
                                                                        indel_realign_job.rv(6))
    
    # Run analysis workflow
    if options.run_analysis:
        joined_sibling_names = [proband_name]
        if siblings_names is not None: 
            joined_sibling_names += siblings_names
        analysis_workflow_job = stage4_jobs.addFollowOnJobFn(run_analysis, context, final_pedigree_joint_called_vcf,
                                                           final_maternal_bam, final_maternal_bam_index, 
                                                           final_paternal_bam, final_paternal_bam_index,
                                                           final_sibling_bam_list_job.rv(0), final_sibling_bam_list_job.rv(1),
                                                           proband_name, maternal_name, paternal_name,
                                                           joined_sibling_names, options.sibling_genders, options.sibling_affected,
                                                           options.bypass, options.cadd_lines,
                                                           options.chrom_dir, options.edit_dir,
                                                           options.split_lines, options.genome_build, options.cadd_data)
        return final_pedigree_joint_called_vcf, final_pedigree_joint_called_vcf_index, final_proband_bam, final_proband_bam_index, final_proband_gvcf, final_proband_gvcf_index, final_maternal_bam, final_maternal_bam_index, final_maternal_gvcf, final_maternal_gvcf_index, final_paternal_bam, final_paternal_bam_index, final_paternal_gvcf, final_paternal_gvcf_index, final_sibling_bam_list_job.rv(0), final_sibling_bam_list_job.rv(1), final_sibling_gvcf_list, final_sibling_gvcf_index_list, trio_joint_called_vcf, trio_joint_called_vcf_index, analysis_workflow_job.rv()
    else:
        return final_pedigree_joint_called_vcf, final_pedigree_joint_called_vcf_index, final_proband_bam, final_proband_bam_index, final_proband_gvcf, final_proband_gvcf_index, final_maternal_bam, final_maternal_bam_index, final_maternal_gvcf, final_maternal_gvcf_index, final_paternal_bam, final_paternal_bam_index, final_paternal_gvcf, final_paternal_gvcf_index, final_sibling_bam_list_job.rv(0), final_sibling_bam_list_job.rv(1), final_sibling_gvcf_list, final_sibling_gvcf_index_list, trio_joint_called_vcf, trio_joint_called_vcf_index

def pedigree_main(context, options):
    """
    Wrapper for vg pedigree. 
    """

    validate_pedigree_options(context, options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Make an index collection
            indexes = {}
           
            # Upload each index we have
            if options.xg_index is not None:
                indexes['xg'] = importer.load(options.xg_index)
            if options.gcsa_index is not None:
                indexes['gcsa'] = importer.load(options.gcsa_index)
                indexes['lcp'] = importer.load(options.gcsa_index + ".lcp")
            if options.gbwt_index is not None:
                indexes['gbwt'] = importer.load(options.gbwt_index)
            if options.distance_index is not None:
                indexes['distance'] = importer.load(options.distance_index)
            if options.minimizer_index is not None:
                indexes['minimizer'] = importer.load(options.minimizer_index)
            if options.snarls_index is not None:
                indexes['snarls'] = importer.load(options.snarls_index)
            if options.id_ranges is not None:
                indexes['id_ranges'] = importer.load(options.id_ranges)
            
            # Upload ref fasta files to the remote IO Store
            # ref_fasta_id, ref_fasta_index_id, ref_fasta_dict_id
            inputRefFastaFileIDs = []
            inputRefFastaFileIDs.append(importer.load(options.ref_fasta))
            inputRefFastaFileIDs.append(importer.load(options.ref_fasta_index))
            inputRefFastaFileIDs.append(importer.load(options.ref_fasta_dict))
            
            # Upload miscelaneous file required by pedigree workflow to the remote IO Store
            # path_list_id, ped_file, snpeff_annotation_database, genetic_map
            inputMiscFileIDs = {}
            if options.path_list is not None:
                inputMiscFileIDs['path_list'] = importer.load(options.path_list)
            if options.ped_file is not None:
                inputMiscFileIDs['ped_file'] = importer.load(options.ped_file)
            if options.snpeff_database is not None:
                inputMiscFileIDs['snpeff_data'] = importer.load(options.snpeff_database)
            if options.genetic_map is not None:
                inputMiscFileIDs['genetic_map'] = importer.load(options.genetic_map)
            
            # Upload other local files to the remote IO Store
            inputReadsFileIDsProband = []
            if options.fastq_proband:
                for sample_reads in options.fastq_proband:
                    inputReadsFileIDsProband.append(importer.load(sample_reads))
            elif options.gam_input_reads_proband:
                inputReadsFileIDsProband.append(importer.load(options.gam_input_reads_proband))
            else:
                assert options.bam_input_reads_proband
                inputReadsFileIDsProband.append(importer.load(options.bam_input_reads_proband))

            inputReadsFileIDsMaternal = []
            if options.fastq_maternal:
                for sample_reads in options.fastq_maternal:
                    inputReadsFileIDsMaternal.append(importer.load(sample_reads))
            elif options.gam_input_reads_maternal:
                inputReadsFileIDsMaternal.append(importer.load(options.gam_input_reads_maternal))
            else:
                assert options.bam_input_reads_maternal
                inputReadsFileIDsMaternal.append(importer.load(options.bam_input_reads_maternal))
            
            inputReadsFileIDsPaternal = []
            if options.fastq_paternal:
                for sample_reads in options.fastq_paternal:
                    inputReadsFileIDsPaternal.append(importer.load(sample_reads))
            elif options.gam_input_reads_paternal:
                inputReadsFileIDsPaternal.append(importer.load(options.gam_input_reads_paternal))
            else:
                assert options.bam_input_reads_paternal
                inputReadsFileIDsPaternal.append(importer.load(options.bam_input_reads_paternal))
            
            inputReadsFileIDsSiblings = []
            if options.fastq_siblings:
                for sample_reads in options.fastq_siblings:
                    inputReadsFileIDsSiblings.append(importer.load(sample_reads))
            elif options.gam_input_reads_siblings:
                for sample_gam_reads in options.gam_input_reads_siblings:
                    inputReadsFileIDsSiblings.append(importer.load(sample_gam_reads))
            elif options.bam_input_reads_siblings:
                for sample_bam_reads in options.bam_input_reads_siblings:
                    inputReadsFileIDsSiblings.append(importer.load(sample_bam_reads))
            
            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_pedigree, context, options, options.fastq_proband,
                                     options.gam_input_reads_proband, options.bam_input_reads_proband,
                                     options.fastq_maternal,
                                     options.gam_input_reads_maternal, options.bam_input_reads_maternal,
                                     options.fastq_paternal,
                                     options.gam_input_reads_paternal, options.bam_input_reads_paternal,
                                     options.fastq_siblings,
                                     options.gam_input_reads_siblings, options.bam_input_reads_siblings,
                                     options.proband_name,
                                     options.maternal_name,
                                     options.paternal_name,
                                     options.sibling_names,
                                     options.indel_realign_bams, options.snpeff_annotation,
                                     options.interleaved, options.mapper,
                                     importer.resolve(indexes),
                                     importer.resolve(inputRefFastaFileIDs),
                                     importer.resolve(inputMiscFileIDs),
                                     reads_file_ids_proband=importer.resolve(inputReadsFileIDsProband),
                                     reads_file_ids_maternal=importer.resolve(inputReadsFileIDsMaternal),
                                     reads_file_ids_paternal=importer.resolve(inputReadsFileIDsPaternal),
                                     reads_file_ids_siblings=importer.resolve(inputReadsFileIDsSiblings),
                                     bam_output=options.bam_output, surject=options.surject,
                                     validate=options.validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

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
    
