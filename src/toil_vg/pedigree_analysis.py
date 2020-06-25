#!/usr/bin/env python
"""
pedigree_analysis.py: run the UDP analysis workflow on a pedigree processed by the 'toil-vg pedigree' subcommand workflow.

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

def analysis_subparser(parser):
    """
    Create a subparser for the analysis workflow.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    # VCFtoShebang options
    parser.add_argument("--cohort_vcf", type=make_url, required=True,
        help="Path to cohort vcf processed by toil-vg pedigree")    
    parser.add_argument("--sample_name", type=str, required=True,
        help="Cohort proband sample name as present in the sample column of --cohort_vcf")
    parser.add_argument("--bypass", action="store_true", default=False,
        help="Parameter for vcftoshebang.")
    # TODO: what is the bypass parameter?
    parser.add_argument("--cadd_lines", type=int, default=4985,
        help="Parameter for vcftoshebang.")
    # TODO: what is the cadd lines parameter?
    parser.add_argument("--chrom_dir", type=str, required=True,
        help="Path to chromosome annotation directory used by vcftoshebang")
    parser.add_argument("--edit_dir", type=str, required=True,
        help="Path to directory containing master edit files used by vcftoshebang")
    
    # CADD options
    parser.add_argument("--split_lines", type=int, default=30000,
        help="Number of lines to chunk the input VCF for CADD processing.")
    parser.add_argument("--genome_build", type=str, default="GRCh37",
        help="Genome annotation version for the CADD engine")
    parser.add_argument("--cadd_data", type=str, required=True,
        help="Path to cadd engine data directory")
    
    # BMTB options
    parser.add_argument("--maternal_bam", type=make_url, required=True,
        help="Input maternal bam file")    
    parser.add_argument("--maternal_bai", type=make_url, required=True,
        help="Input maternal bam file .bai index")    
    parser.add_argument("--paternal_bam", type=make_url, required=True,
        help="Input paternal bam file")    
    parser.add_argument("--paternal_bai", type=make_url, required=True,
        help="Input paternal bam file .bai index")    
    parser.add_argument("--siblings_bam", nargs='+', type=make_url, required=True,
        help="INPUT sibling bam files. Proband bam file must be listed 1st.")
    parser.add_argument("--siblings_bai", nargs='+', type=make_url, required=True,
        help="INPUT sibling bam file .bai indicies. Same sample order as --siblings_bam.")
    parser.add_argument("--maternal_name", type=str, required=True,
        help="Sample name of the pedigree mother.")
    parser.add_argument("--paternal_name", type=str, required=True,
        help="Sample name of the pedigree father.")
    parser.add_argument("--sibling_names", nargs='+', type=str, required=True,
        help="Sample names of the siblings. Same sample order as --siblings_bam.")
    parser.add_argument("--sibling_genders", nargs='+', type=int, required=True,
        help="Gender of each sibling sample. 0 = male, 1 = female. Same sample order as --siblings_bam.")
    parser.add_argument("--sibling_affected", nargs='+', type=int, required=True,
        help="Affected status of each sibling sample. 0 = unaffected, 1 = affected. Same sample order as --siblings_bam.")
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def run_vcftoshebang(job, context, sample_name, cohort_vcf_id, bypass, cadd_lines, chrom_dir, edit_dir):
    """ run vcftoshebang on a cohort vcf.
    
    Takes a joint-called cohort vcf from the vg_pedigree.py workflow along with the sample name of the proband,
    the bypass boolean parameter for vcftoshebang, the number of lines to split the varsifter variant list that
    will be passed to the CADD engine, the full path to the directory containing chromosome config files, and
    the full path to the directory containing the master edit files.
    
    Returns a tuple containing varsifter formatted and annotated version of the cohort vcf as the first element, and
    a list of file ids that will passed to the CADD engine if any exist.
    
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    cohort_vcf_file = os.path.join(work_dir, os.path.basename(cohort_vcf_id))
    job.fileStore.readGlobalFile(cohort_vcf_id, cohort_vcf_file)
    
    # Decompress cohort vcf if already compressed
    if os.path.basename(cohort_vcf_file).endswith('.gz'):
        context.runner.call(job, ['bgzip', '-d', os.path.basename(cohort_vcf_file)], work_dir = work_dir, tool_name='vg')
        cohort_vcf_file = os.path.splitext(cohort_vcf_file)[0]
    
    # Remove the hs37d5 contig
    context.runner.call(job, ['vcftools', '--vcf', os.path.basename(cohort_vcf_file), '--not-chr', 'hs37d5', '--recode-INFO-all', '--recode', '--out', '{}.filtered'.format(os.path.basename(os.path.splitext(cohort_vcf_file)[0]))], work_dir = work_dir, tool_name='vcftools')
    input_vcf_file = "{}.filtered.recode.vcf".format(os.path.basename(os.path.splitext(cohort_vcf_file)[0]))
    
    output_dir = "vcf2shebang_output/"
    bypass_conf = "NO"
    if bypass == 'true': bypass_conf = "YES"
    context.runner.call(job, ['mkdir', '{}'.format(output_dir)], work_dir = work_dir)
    cmd_list = []
    cmd_list.append(['cp', '/vcftoshebang/VCFtoShebang_Config.txt', '.'])
    cmd_list.append(['sed', '-i', '\"s|.*PROBAND_NAME.*|PROBAND_NAME\t{}|\"'.format(sample_name), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*OUTPUT_DIR.*|OUTPUT_DIR\t{}|\"'.format(output_dir), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*UNROLLED_VCF_PATH.*|UNROLLED_VCF_PATH\t{}|\"'.format(input_vcf_file), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*BYPASS.*|BYPASS\t{}|\"'.format(bypass_conf), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*CADD_LINES.*|CADD_LINES\t{}|\"'.format(cadd_lines), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*CHROM_DIR.*|CHROM_DIR\t$PWD/{}|\"'.format(os.path.basename(os.path.normpath(chrom_dir))), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*EDIT_DIR.*|EDIT_DIR\t$PWD/{}|\"'.format(os.path.basename(os.path.normpath(edit_dir))), 'VCFtoShebang_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*EDITOR_CONFIG.*|EDITOR_CONFIG\t/vcftoshebang/edit_config.txt|\"', 'VCFtoShebang_Config.txt'])
    cmd_list.append(['java', '-XX:+UnlockExperimentalVMOptions', '-XX:ActiveProcessorCount=32', '-cp', '/vcftoshebang/VCFtoShebang.jar:/vcftoshebang/json_simple.jar',
                             'Runner', 'VCFtoShebang_Config.txt'])
    chain_cmds = [' '.join(p) for p in cmd_list]
    command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' && '.join(chain_cmds))]
    context.runner.call(job, command, work_dir = work_dir, tool_name='vcf2shebang', mount_list=[chrom_dir,edit_dir])
    is_empty = True
    with open(os.path.join(work_dir, 'vcf2shebang_output/{}_unrolled_snpeff_fix_overlap_mono_CADD_Input_Files/{}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt'.format(sample_name,sample_name)), 'r') as cadd_input_file:
        for line in cadd_input_file:
            if "#" in line: continue
            else:
                is_empty = False
                break
    
    output_vs_path = os.path.join(work_dir, 'vcf2shebang_output/{}_unrolled_snpeff_fix_overlap_mono_shebang.vs'.format(sample_name))
    output_cadd_vcf_path = None
    if not is_empty:
        output_cadd_vcf_path = os.path.join(work_dir, 'vcf2shebang_output/{}_unrolled_snpeff_fix_overlap_mono_CADD_Input_Files/{}_unrolled_snpeff_fix_overlap_mono_CADD_input_file.txt.gz'.format(sample_name,sample_name))
    
    return context.write_output_file(job, output_vs_path), context.write_output_file(job, output_cadd_vcf_path)
    
def run_split_vcf(job, context, vcf_file_id, split_lines):
    """ split vcf into chunks for passing through to the CADD engine for CADD scoring.
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_file_id))
    job.fileStore.readGlobalFile(vcf_file_id, vcf_file)
    
    tempname = "{}_tmp".format(os.path.basename(os.path.splitext(os.path.splitext(vcf_file)[0])[0]))
    cmd = [['zcat', os.path.basename(vcf_file)]]
    cmd.append(['grep', '-v', '^GL|^#'])
    cmd.append(['cut', '-f', '1-5'])
    cmd.append(['split', '-l', str(split_lines), '--additional-suffix=.vcf', '-d', '-', tempname])
    context.runner.call(job, cmd, work_dir = work_dir)
    
    vcf_chunk_ids = []
    for chunk_name in sorted(os.listdir(work_dir)):
        if chunk_name.endswith('.vcf') and 'tmp' in chunk_name:
            vcf_chunk_ids.append(context.write_intermediate_file(job, os.path.join(work_dir, chunk_name)))
    
    return vcf_chunk_ids

def run_cadd(job, context, chunk_vcf_id, genome_build, cadd_data_dir):
    """ run the main CADD engine on vcf chunks.
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(chunk_vcf_id))
    job.fileStore.readGlobalFile(chunk_vcf_id, vcf_file)
    
    base_vcf_name = os.path.basename(os.path.splitext(vcf_file)[0])
    cadd_data_dir_basename = os.path.basename(cadd_data_dir)
    cmd_list = []
    cmd_list.append(['source', 'activate', '$(head', '-1', '/usr/src/app/environment.yml', '|', 'cut', '-d\'', '\'', '-f2)'])
    cmd_list.append(['/bin/bash', '/usr/src/app/CADD.sh', '-v', '\"v1.5\"', '-g', genome_build, '-o', '$PWD/{}_out.tsv.gz'.format(base_vcf_name), '-d', '$PWD/{}'.format(cadd_data_dir_basename), os.path.basename(vcf_file)])
    chain_cmds = [' '.join(p) for p in cmd_list]
    command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' && '.join(chain_cmds))]
    context.runner.call(job, command, work_dir = work_dir, tool_name='cadd', mount_list=[cadd_data_dir])
    
    output_cadd_path = os.path.join(work_dir, '{}_out.tsv.gz'.format(base_vcf_name))
    return context.write_output_file(job, output_cadd_path)

def run_merge_annotated_vcf(job, context, cadd_output_chunk_ids):
    """ run the merge task for the CADD engine output chunks.
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    cmd_list = [['source', 'activate', '$(head', '-1', '/usr/src/app/environment.yml', '|', 'cut', '-d\' \'', '-f2)']]
    # Decompress and concatenate cadd output chunks into a single file
    cadd_chunk_merged_file = os.path.join(work_dir, 'merged_CADDv1.5_offline_unsorted')
    for cadd_output_chunk_id in sorted(cadd_output_chunk_ids):
        cadd_chunk_file = os.path.join(work_dir, os.path.basename(cadd_output_chunk_id))
        job.fileStore.readGlobalFile(cadd_output_chunk_id, cadd_chunk_file)
        cmd_list.append(['zcat', os.path.basename(cadd_chunk_file), '>>', 'merged_CADDv1.5_offline_unsorted'])
        
    # Sort and process the merged file
    cmd_list.append(['sort', '-k1,1', '-k2,2n', 'merged_CADDv1.5_offline_unsorted', '>', 'merged_CADDv1.5_offline.vcf'])
    cmd_list.append(['rm', '-f', 'merged_CADDv1.5_offline_unsorted'])
    cmd_list.append(['python', '/usr/src/app/CADD_offline_mito_postprocessing.py', '-c', '/usr/src/app/whole_mito_SNP_pp2_predictions_sorted.txt', '-i', 'merged_CADDv1.5_offline.vcf', '-o', 'merged_CADDv1.5_offline_proper_format.vcf'])
    cmd_list.append(['rm', '-f', 'merged_CADDv1.5_offline.vcf'])
    chain_cmds = [' '.join(p) for p in cmd_list]
    command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' && '.join(chain_cmds))]
    context.runner.call(job, command, work_dir = work_dir, tool_name='cadd')
    
    merged_cadd_output_vcf_path = os.path.join(work_dir, 'merged_CADDv1.5_offline_proper_format.vcf')
    return context.write_output_file(job, merged_cadd_output_vcf_path)

def run_cadd_editor(job, context, vcftoshebang_vs_file_id, merged_cadd_vcf_file_id):
    """ run editor for merged CADD output postprocessing and merges annotations with vcftoshebang varsifter output.
        Outputs the file in Varsifter format.
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vcftoshebang_vs_file = os.path.join(work_dir, os.path.basename(vcftoshebang_vs_file_id))
    job.fileStore.readGlobalFile(vcftoshebang_vs_file_id, vcftoshebang_vs_file)
    
    merged_cadd_vcf_file = os.path.join(work_dir, os.path.basename(merged_cadd_vcf_file_id))
    job.fileStore.readGlobalFile(merged_cadd_vcf_file_id, merged_cadd_vcf_file)
    
    command = ['java', '-cp', '/cadd_edit/NewCaddEditor.jar:/cadd_edit/commons-cli-1.4.jar', 'NewCaddEditor',
                '--input_vs', os.path.basename(vcftoshebang_vs_file), '--output_cadd', os.path.basename(merged_cadd_vcf_file), '--output_vs', 'cadd_editor_output.vs']
    context.runner.call(job, command, work_dir = work_dir, tool_name='caddeditor')
    
    cadd_editor_output_path = os.path.join(work_dir, 'cadd_editor_output.vs')
    return context.write_output_file(job, cadd_editor_output_path)

def run_bmtb(job, context, analysis_ready_vs_file_id,
               maternal_bam_id, maternal_bai_id, paternal_bam_id, paternal_bai_id, sibling_bam_ids, sibling_bai_ids,
               maternal_name, paternal_name, sibling_names, sibling_genders, sibling_affected):
    """ run the Black Magic Toolbox candidate analysis program on the input varsifter file, cohort bam files, and sibling gender
        and sibling affected status.
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vs_file_path = os.path.join(work_dir, os.path.basename(analysis_ready_vs_file_id))
    job.fileStore.readGlobalFile(analysis_ready_vs_file_id, vs_file_path)
    
    m_bam_path = os.path.join(work_dir, os.path.basename(maternal_bam_id))
    job.fileStore.readGlobalFile(maternal_bam_id, m_bam_path)
    m_bai_path = os.path.join(work_dir, os.path.basename(maternal_bai_id))
    job.fileStore.readGlobalFile(maternal_bai_id, m_bai_path)
    
    f_bam_path = os.path.join(work_dir, os.path.basename(paternal_bam_id))
    job.fileStore.readGlobalFile(paternal_bam_id, f_bam_path)
    f_bai_path = os.path.join(work_dir, os.path.basename(paternal_bai_id))
    job.fileStore.readGlobalFile(paternal_bai_id, f_bai_path)
    
    s_bam_paths = []
    for s_bam_id in sibling_bam_ids:
        s_bam_path = os.path.join(work_dir, os.path.basename(s_bam_id))
        job.fileStore.readGlobalFile(s_bam_id, s_bam_path)
        s_bam_paths.append(s_bam_path)
    s_bai_paths = []
    for s_bai_id in sibling_bai_ids:
        s_bai_path = os.path.join(work_dir, os.path.basename(s_bai_id))
        job.fileStore.readGlobalFile(s_bai_id, s_bai_path)
        s_bai_paths.append(s_bai_path)
    
    context.runner.call(job, ['ls', '-l', '/bmtb/'], work_dir = work_dir, tool_name='bmtb')
    cmd_list = [['cp', '-r', '/bmtb/Configs', '$PWD/Configs']]
    cmd_list.append(['rm', '-f', '$PWD/Configs/BAM_Directory_Config.txt'])
    cmd_list.append(['touch', '$PWD/Configs/BAM_Directory_Config.txt'])
    cmd_list.append(['echo', '-e', '\"{}\t$PWD/{}\"'.format(maternal_name,os.path.basename(m_bam_path)), '>>', '$PWD/Configs/BAM_Directory_Config.txt'])
    cmd_list.append(['echo', '-e', '\"{}\t$PWD/{}\"'.format(paternal_name,os.path.basename(f_bam_path)), '>>', '$PWD/Configs/BAM_Directory_Config.txt'])
    sibling_id_string = "NA"
    for i, (s_bam_path,s_bai_path,s_name,s_gender,s_affected) in enumerate(zip(s_bam_paths,s_bai_paths,sibling_names,sibling_genders,sibling_affected)):
        cmd_list.append(['echo', '-e', '\"{}\t$PWD/{}\"'.format(s_name,os.path.basename(s_bam_path)), '>>', '$PWD/Configs/BAM_Directory_Config.txt'])
        # Add proband data to master config file
        if i == 0:
            cmd_list.append(['sed', '-i', '\"s|.*PB_ID.*|PB_ID\t{}|\"'.format(s_name), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
            cmd_list.append(['sed', '-i', '\"s|.*PB_GENDER.*|PB_GENDER\t{}|\"'.format(s_gender), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
            cmd_list.append(['sed', '-i', '\"s|.*OUT_FILE_NAME.*|OUT_FILE_NAME\t{}|\"'.format(s_name), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
        else:
            if sibling_id_string == "NA":
                sibling_id_string = "{},{},{};".format(s_name,s_gender,s_affected)
            else:
                sibling_id_string += "{},{},{};".format(s_name,s_gender,s_affected)
    
    cmd_list.append(['sed', '-i', '\"s|.*VS_FILE_PATH.*|VS_FILE_PATH\t$PWD/{}|\"'.format(os.path.basename(vs_file_path)), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*CONFIG_FILE_DIREC.*|CONFIG_FILE_DIREC\t$PWD/Configs|\"', '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*FATHER_ID.*|FATHER_ID\t{}|\"'.format(paternal_name), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*MOTHER_ID.*|MOTHER_ID\t{}|\"'.format(maternal_name), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
    cmd_list.append(['sed', '-i', '\"s|.*SIB_IDS.*|SIB_IDS\t{}|\"'.format(sibling_id_string), '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
    cmd_list.append(['java', '-cp', '/bmtb/bmtb.jar:/bmtb/htsjdk-2.19.0-47-gc5ed6b7-SNAPSHOT.jar', 'general.Runner', '$PWD/Configs/BMTB_Genome_Input_Config.txt'])
    cmd_list.append(['tar', 'czvf', '\"{}_BlackBox_Output.tar.gz\"'.format(sibling_names[0]), '\"{}_BlackBox_Output\"'.format(sibling_names[0])])
    chain_cmds = [' '.join(p) for p in cmd_list]
    command = ['/bin/bash', '-c', 'set -eo pipefail && {}'.format(' && '.join(chain_cmds))]
    context.runner.call(job, command, work_dir = work_dir, tool_name='bmtb')
    
    output_package_path = os.path.join(work_dir, '{}_BlackBox_Output.tar.gz'.format(sibling_names[0]))
    return context.write_output_file(job, output_package_path)

def run_cadd_jobs(job, context, vcf_chunk_ids, genome_build, cadd_data_dir):
    """ helper function for running multiple cadd jobs per vcf chunk
    """
    cadd_engine_output_ids = []
    for vcf_chunk_id in vcf_chunk_ids:
        cadd_job = job.addChildJobFn(run_cadd, context, vcf_chunk_id, genome_build, cadd_data_dir,
                                        cores=context.config.misc_cores,
                                        memory=context.config.misc_mem,
                                        disk=context.config.misc_disk)
        cadd_engine_output_ids.append(cadd_job.rv())
    
    return cadd_engine_output_ids
    
def run_analysis(job, context, cohort_vcf_id,
                       maternal_bam_id, maternal_bai_id, paternal_bam_id, paternal_bai_id, sibling_bam_ids, sibling_bai_ids,
                       sample_name, maternal_name, paternal_name,
                       sibling_names, sibling_genders, sibling_affected,
                       bypass, cadd_lines,
                       chrom_dir, edit_dir,
                       split_lines, genome_build, cadd_data_dir):
    """ run vcf to shebang varsifter file conversion, then do cadd scoring and annotation, finally run the blackmagiktoolbox workflow.
        returns final candidate varsifter file, paired with total surject time
    (excluding toil-vg overhead such as transferring and splitting files )"""

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)
    
    vcf_to_shebang_job = child_job.addChildJobFn(run_vcftoshebang, context, sample_name, cohort_vcf_id, bypass, cadd_lines, chrom_dir, edit_dir,
                                                    cores=context.config.alignment_cores,
                                                    memory=context.config.alignment_mem,
                                                    disk=context.config.alignment_disk)
    
    analysis_ready_vs_file_id = vcf_to_shebang_job.rv(0)
    if vcf_to_shebang_job.rv(1) is not None:
        RealtimeLogger.info("Some variants don't have CADD scores, running them through the CADD engine workflow.")
        split_vcf_job = vcf_to_shebang_job.addChildJobFn(run_split_vcf, context, vcf_to_shebang_job.rv(1), cadd_lines)
        cadd_jobs = split_vcf_job.addFollowOnJobFn(run_cadd_jobs, context, split_vcf_job.rv(), genome_build, cadd_data_dir,
                                                cores=context.config.misc_cores,
                                                memory=context.config.misc_mem,
                                                disk=context.config.misc_disk)
        merge_annotated_vcf_job = cadd_jobs.addFollowOnJobFn(run_merge_annotated_vcf, context, cadd_jobs.rv())
        cadd_edit_job = merge_annotated_vcf_job.addFollowOnJobFn(run_cadd_editor, context, vcf_to_shebang_job.rv(0), merge_annotated_vcf_job.rv(),
                                                                    cores=context.config.misc_cores,
                                                                    memory=context.config.alignment_mem,
                                                                    disk=context.config.alignment_disk)
        
        analysis_ready_vs_file_id = cadd_edit_job.rv()
    
    bmtb_job = vcf_to_shebang_job.addFollowOnJobFn(run_bmtb, context, analysis_ready_vs_file_id,
                                                   maternal_bam_id, maternal_bai_id, paternal_bam_id, paternal_bai_id, sibling_bam_ids, sibling_bai_ids,
                                                   maternal_name, paternal_name, sibling_names, sibling_genders, sibling_affected,
                                                   cores=context.config.misc_cores,
                                                   memory=context.config.alignment_mem,
                                                   disk=context.config.alignment_disk)
    
    return bmtb_job.rv()
    
def analysis_main(context, options):
    """
    Wrapper for the pedigree analysis pipeline. 
    """

    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Upload local files to the remote IO Store
            inputVCFFileID = importer.load(options.cohort_vcf)
            inputMBAMFileID = importer.load(options.maternal_bam)
            inputMBAMINDEXFileID = importer.load(options.maternal_bai)
            inputFBAMFileID = importer.load(options.paternal_bam)
            inputFBAMINDEXFileID = importer.load(options.paternal_bai)
            inputSiblingBAMFileIDs = []
            for sibling_bam in options.siblings_bam:
                inputSiblingBAMFileIDs.append(importer.load(sibling_bam))
            inputSiblingBAMINDEXFileIDs = []
            for sibling_bai in options.siblings_bai:
                inputSiblingBAMINDEXFileIDs.append(importer.load(sibling_bai))
            
            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_analysis, context,
                                     importer.resolve(inputVCFFileID),
                                     importer.resolve(inputMBAMFileID),
                                     importer.resolve(inputMBAMINDEXFileID),
                                     importer.resolve(inputFBAMFileID),
                                     importer.resolve(inputFBAMINDEXFileID),
                                     importer.resolve(inputSiblingBAMFileIDs),
                                     importer.resolve(inputSiblingBAMINDEXFileIDs),
                                     options.sample_name,
                                     options.maternal_name,
                                     options.paternal_name,
                                     options.sibling_names, options.sibling_genders, options.sibling_affected,
                                     options.bypass, options.cadd_lines,
                                     options.chrom_dir, options.edit_dir,
                                     options.split_lines, options.genome_build, options.cadd_data,
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
    
