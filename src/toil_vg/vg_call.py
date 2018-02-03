#!/usr/bin/env python2.7
"""
Generate a VCF from a GAM and XG by splitting into GAM/VG chunks.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno, copy, time
from uuid import uuid4
import logging
from collections import defaultdict

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore

logger = logging.getLogger(__name__)

def call_subparser(parser):
    """
    Create a subparser for calling.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("xg_path", type=make_url,
                        help="input xg file")
    parser.add_argument("sample_name", type=str,
                        help="sample name (ex NA12878)")
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified "
                        "using same syntax as toil jobStore")
    parser.add_argument("--chroms", nargs='+', required=True,
                        help="Name(s) of reference path in graph(s) (separated by space)."
                        " Must be same length/order as --gams")
    # todo: move to chunked_call_parse_args and share with toil-vg run
    parser.add_argument("--gams", nargs='+', required=True, type=make_url,
                        help="GAMs to call.  One per chromosome. Must be same length/order as --chroms")
    parser.add_argument("--gam_index_cores", type=int,
                        help="number of threads used for gam indexing")
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common calling options shared with toil_vg pipeline
    chunked_call_parse_args(parser)

    # Add common docker options shared with toil_vg pipeline
    add_container_tool_parse_args(parser)
                        

def chunked_call_parse_args(parser):
    """ centralize calling parameters here """
    parser.add_argument("--overlap", type=int,
                        help="overlap option that is passed into make_chunks and call_chunk")
    parser.add_argument("--call_chunk_size", type=int,
                        help="chunk size")
    parser.add_argument("--call_opts", type=str,
                        help="arguments to pass to vg call (wrapped in \"\")")
    parser.add_argument("--genotype", action="store_true",
                        help="use vg genotype instead of vg call")
    parser.add_argument("--no_augment", action="store_true",
                        help="only call existing variants in the graph")
    parser.add_argument("--genotype_opts", type=str,
                        help="arguments to pass to vg genotype (wrapped in \"\")")
    parser.add_argument("--filter_opts", type=str,
                        help="argument to pass to vg filter (wrapped in \"\")")
    parser.add_argument("--calling_cores", type=int,
                        help="number of threads during the variant calling step")
    parser.add_argument("--call_chunk_cores", type=int,
                        help="number of threads used for extracting chunks for calling")
    parser.add_argument("--vcf_offsets", nargs='+', default=[],
                         help="offset(s) to apply to output vcfs(s). (order of --chroms)")

def validate_call_options(options):    
    require(len(options.chroms) == len(options.gams) or len(options.gams) == 1,
            'Number of --chroms must be 1 or same as number of --gams')
    require(not options.vcf_offsets or len(options.vcf_offsets) == len(options.chroms),
            'Number of --vcf_offsets if specified must be same as number of --chroms')
    require(not options.no_augment or options.genotype,
            '--no_augment currently only supported with --genotype')
    
def sort_vcf(job, drunner, vcf_path, sorted_vcf_path):
    """ from vcflib """
    vcf_dir, vcf_name = os.path.split(vcf_path)
    with open(sorted_vcf_path, "w") as outfile:
        drunner.call(job, [['bcftools', 'view', '-h', vcf_name]], outfile=outfile,
                     work_dir=vcf_dir)
    with open(sorted_vcf_path, "a") as outfile:
        drunner.call(job, [['bcftools', 'view', '-H', vcf_name],
                      ['sort', '-k1,1d', '-k2,2n']], outfile=outfile,
                     work_dir=vcf_dir)

def run_vg_call(job, context, sample_name, vg_id, gam_id, xg_id = None,
                path_names = [], seq_names = [], seq_offsets = [], seq_lengths = [],
                filter_opts = [], augment_opts = [], call_opts = [],
                keep_pileup = False, keep_xg = False, keep_gam = False,
                keep_augmented = False, chunk_name = 'call', genotype = False,
                augment = True):
    """ Run vg call on a single graph.

    NOTE: Can now run vg genotype as well, but the plan is to fold genotype into call
    soon.  At that point, will remove run_vg_genotype() function entirely, which is kinda
    sitting around useless now. 

    Returns (vcf_id, pileup_id, xg_id, gam_id, augmented_graph_id).  pileup_id and xg_id
    can be same as input if they are not computed.  If pileup/xg/augmented are 
    computed, the returned ids will be None unless appropriate keep_flag set
    (to prevent sending them to the file store if they aren't wanted)

    User is responsible to make sure that options passed in context.config.*_opts don't conflict
    with seq_names, seq_offsets, seq_lengths etc. If not provided, the pileup is computed.

    gam filtering is only done if filter_opts are passed in. 

    chunk_name option is only for working filenames (to make more readable)

    """

    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    vg_path = os.path.join(work_dir, '{}.vg'.format(chunk_name))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(chunk_name))
    job.fileStore.readGlobalFile(gam_id, gam_path)
    xg_path = os.path.join(work_dir, '{}.xg'.format(chunk_name))
    defray = filter_opts and ('-D' in filter_opts or '--defray-ends' in filter_opts)
    if xg_id and defray:
        job.fileStore.readGlobalFile(xg_id, xg_path)
    pu_path = os.path.join(work_dir, '{}.pu'.format(chunk_name))
    trans_path = os.path.join(work_dir, '{}.trans'.format(chunk_name))
    support_path = os.path.join(work_dir, '{}.support'.format(chunk_name))
    aug_path = os.path.join(work_dir, '{}_aug.vg'.format(chunk_name))
    aug_gam_path = os.path.join(work_dir, '{}_aug.gam'.format(chunk_name))

    timer = TimeTracker()

    # we only need an xg if using vg filter -D
    if not xg_id and defray:
        timer.start('chunk-xg')
        context.runner.call(job, ['vg', 'index', os.path.basename(vg_path), '-x',
                                   os.path.basename(xg_path), '-t', str(context.config.calling_cores)],
                             work_dir = work_dir)
        timer.stop()
        if keep_xg:
            xg_id = context.write_intermediate_file(job, xg_path)

    # optional gam filtering
    gam_filter_path = gam_path + '.filter'
    filter_command = None
    if filter_opts:
        filter_command = ['vg', 'filter', os.path.basename(gam_path), '-t', '1'] + filter_opts
        if defray:
            filter_command += ['-x', os.path.basename(xg_path)]

    # we filter separated when running genotype (due to augment -A)
    if filter_command and genotype:
        with open(gam_filter_path, 'w') as gam_filter_stream:
            timer.start('call-filter')
            context.runner.call(job, filter_command, work_dir=work_dir, outfile=gam_filter_stream)
            timer.stop()
        gam_path = gam_filter_path
        filter_command = None
        
    # augment command with optional filter piped at beginning
    augment_generated_opts = ['-Z', os.path.basename(trans_path)]
    if keep_pileup:
        augment_generated_opts += ['-P', os.path.basename(pu_path)]
    if genotype:
        # Make sure to use the augmentation mode that we need for genotype (augment with everything but no supports)
        augment_generated_opts += ['-a', 'direct']
        # We need to keep the augmented gam
        augment_generated_opts += ['-A', os.path.basename(aug_gam_path)]
    else:
        # Make sure to use the augmentation mode for vg call (which can calculate supports)
        augment_generated_opts += ['-a', 'pileup']
        # And calculate the supports instead of the augmented gam
        augment_generated_opts += ['-S', os.path.basename(support_path)]
    augment_command = []
    if filter_command is not None:
        aug_gam_input = '-'
        augment_command.append(filter_command)
        if keep_gam:
            augment_command.append(['tee', os.path.basename(gam_filter_path)])
    else:
        aug_gam_input = os.path.basename(gam_path)
    augment_command.append(['vg', 'augment', os.path.basename(vg_path), aug_gam_input,
                    '-t', str(context.config.calling_cores)] + augment_opts + augment_generated_opts)

    vcf_path = os.path.join(work_dir, '{}_call.vcf'.format(chunk_name))

    # call
    try:
        if augment:
            with open(aug_path, 'w') as aug_stream:
                timer.start('call-filter-augment')
                context.runner.call(job, augment_command, work_dir=work_dir, outfile=aug_stream)
                timer.stop()
        else:
            # hack to skip augmentation
            aug_path = vg_path
            aug_gam_path = gam_path
            
        gam_id, pileup_id, aug_graph_id = None, None, None
        if keep_gam and filter_opts:
            gam_id = context.write_intermediate_file(job, gam_filter_path)
        if keep_pileup:
            pileup_id = context.write_intermediate_file(job, pu_path)
        if keep_augmented:
            aug_graph_id = context.write_intermediate_file(job, aug_path)
        
        with open(vcf_path, 'w') as vgcall_stdout:
            if not genotype:
                command = ['vg', 'call', os.path.basename(aug_path), '-t',
                           str(context.config.calling_cores), '-S', sample_name,
                           '-z', os.path.basename(trans_path),
                           '-s', os.path.basename(support_path),
                           '-b', os.path.basename(vg_path)]
            else:
                command = ['vg', 'genotype', os.path.basename(aug_path), '-t',
                           str(context.config.calling_cores), '-s', sample_name,
                           '-v', '-E', '-G', os.path.basename(aug_gam_path)]
            if call_opts:
                command += call_opts
            for path_name in path_names:
                command += ['-r', path_name]
            for seq_name in seq_names:
                command += ['-c', seq_name]
            for seq_length in seq_lengths:
                command += ['-l', seq_length]
            for seq_offset in seq_offsets:
                command += ['-o', seq_offset]

            timer.start('genotype' if genotype else 'call')
            context.runner.call(job, command, work_dir=work_dir,
                                 outfile=vgcall_stdout)
            timer.stop()            

        vcf_id = context.write_intermediate_file(job, vcf_path)

    except Exception as e:
        logging.error("Calling failed. Dumping files.")
        for dump_path in [vg_path, pu_path, gam_filter_path,
                          aug_path, support_path, trans_path, aug_gam_path]:
            if dump_path and os.path.isfile(dump_path):
                context.write_output_file(job, dump_path)        
        raise e
        
    return vcf_id, pileup_id, xg_id, gam_id, aug_graph_id, timer


def run_vg_genotype(job, context, sample_name, vg_id, gam_id, xg_id = None,
                    path_names = [], seq_names = [], seq_offsets = [],
                    seq_lengths = [], genotype_opts = [], filter_opts = [],
                    keep_xg = False, keep_gam = False, keep_augmented = False,
                    index_gam = True, chunk_name = 'genotype'):
    """ Run vg genotype on a single graph.

    Returns (vcf_id, xg_id, gam_id, augmented_graph_id).  xg_id
    can be same as input if not computed.  If xg/augmented are 
    computed, the returned ids will be None unless appropriate keep_flag set
    (to prevent sending them to the file store if they aren't wanted)

    User is responsible to make sure that options passed in context.config.*_opts don't conflict
    with seq_names, seq_offsets, seq_lengths etc. If not provided, the pileup is computed.

    gam filtering is only done if filter_opts are passed in. 

    chunk_name option is only for working filenames (to make more readable)
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    vg_path = os.path.join(work_dir, '{}.vg'.format(chunk_name))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(chunk_name))
    job.fileStore.readGlobalFile(gam_id, gam_path)
    xg_path = os.path.join(work_dir, '{}.xg'.format(chunk_name))
    defray = filter_opts and ('-D' in filter_opts or '--defray-ends' in filter_opts)
    if xg_id and defray:
        job.fileStore.readGlobalFile(xg_id, xg_path)
        
    # we only need an xg if using vg filter -D
    if not xg_id and defray:
        context.runner.call(job, ['vg', 'index', os.path.basename(vg_path), '-x',
                                   os.path.basename(xg_path), '-t', str(context.config.calling_cores)],
                             work_dir = work_dir)
        if keep_xg:
            xg_id = context.write_intermediate_file(job, xg_path)

    # filter
    if filter_opts:
        gam_filter_path = gam_path + '.filter'
        with open(gam_filter_path, 'w') as gam_filter_stream:
            command = [['vg', 'filter', os.path.basename(gam_path), '-t', '1'] + filter_opts]
            if defray:
                command[0] += ['-x', os.path.basename(xg_path)]
            context.runner.call(job, command, work_dir=work_dir, outfile=gam_filter_stream)
        if keep_gam:
            gam_id = context.write_intermediate_file(job, gam_filter_path)
    else:
        gam_filter_path = gam_path    

    # index the gam
    if index_gam:
        # Make a -a index first
        gam_node_index_path = os.path.join(work_dir, gam_filter_path + '.nidx')
        command = ['vg', 'index', '-a', os.path.basename(gam_filter_path),
                   '-d', os.path.basename(gam_node_index_path)]
        context.runner.call(job, command, work_dir=work_dir)

        # Sort the gam
        command = [['vg', 'index', '-A', '-d', os.path.basename(gam_node_index_path)]]        
        # Then make the -N index on the sorted input
        gam_index_path = os.path.join(work_dir, gam_filter_path + '.index')        
        command.append(['vg', 'index', '-N', '-',
                        '-d', os.path.basename(gam_index_path)])
        context.runner.call(job, command, work_dir=work_dir)
    # genotype
    try:
        vcf_path = os.path.join(work_dir, '{}_genotype.vcf'.format(chunk_name))
        aug_graph_id = None
        
        with open(vcf_path, 'w') as vgcall_stdout:
            command = ['vg', 'genotype', os.path.basename(vg_path)]
            if index_gam:
                command += [os.path.basename(gam_index_path)]
            else:
                command += ['-G', os.path.basename(gam_filter_path)]
            command += ['-t', str(context.config.calling_cores)] + genotype_opts + ['-s', sample_name, '-v']
            for path_name in path_names:
                command += ['-r', path_name]
            for seq_name in seq_names:
                command += ['-c', seq_name]
            for seq_length in seq_lengths:
                command += ['-l', seq_length]
            for seq_offset in seq_offsets:
                command += ['-o', seq_offset]
            aug_path = os.path.join(work_dir, '{}_aug.vg'.format(chunk_name))
            if keep_augmented:
                command.append(['-a', os.path.basename(aug_path)])
            context.runner.call(job, command, work_dir=work_dir,
                                 outfile=vgcall_stdout)
            if keep_augmented:
                aug_graph_id = context.write_intermediate_file(job, aug_path)

        vcf_id = context.write_intermediate_file(job, vcf_path)

    except Exception as e:
        logging.error("Failed. Dumping files.")
        context.write_output_file(job, vg_path)
        if not index_gam:
            context.write_output_file(job, gam_filter_path)        
        raise e
        
    return vcf_id, xg_id, gam_id, aug_graph_id
        

def run_call_chunk(job, context, path_name, chunk_i, num_chunks, chunk_offset, clipped_chunk_offset,
                   xg_file_id, vg_chunk_file_id, gam_chunk_file_id, path_size, vcf_offset, sample_name,
                   genotype, augment):
    """ create VCF from a given chunk """

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)

    RealtimeLogger.info("Running call_chunk on path {} and chunk {}".format(path_name, chunk_i))
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Run vg call
    call_job = child_job.addChildJobFn(
        run_vg_call,
        context, sample_name, vg_chunk_file_id, gam_chunk_file_id,
        xg_id = xg_file_id,
        path_names = [path_name], 
        seq_names = [path_name],
        seq_offsets = [chunk_offset + vcf_offset],
        seq_lengths = [path_size],
        filter_opts = context.config.filter_opts,
        augment_opts = context.config.augment_opts,
        call_opts = context.config.call_opts if not genotype else context.config.genotype_opts,
        chunk_name = 'chunk_{}_{}'.format(path_name, chunk_offset),
        genotype = genotype,
        augment = augment,
        cores=context.config.calling_cores,
        memory=context.config.calling_mem, disk=context.config.calling_disk)
    vcf_id, pu_id, xg_id, gam_id, aug_graph_id = [call_job.rv(i) for i in range(5)]
    call_timer = call_job.rv(5)

    clip_job = child_job.addFollowOnJobFn(run_clip_vcf, context, path_name, chunk_i, num_chunks, chunk_offset,
                                          clipped_chunk_offset, vcf_offset, vcf_id,
                                          cores=context.config.calling_cores,
                                          memory=context.config.calling_mem, disk=context.config.calling_disk)

    clip_file_id = clip_job.rv()
    return clip_file_id, call_timer

def run_clip_vcf(job, context, path_name, chunk_i, num_chunks, chunk_offset, clipped_chunk_offset, vcf_offset, vcf_id):
    """ clip the vcf to respect chunk """

     # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # output vcf name
    vcf_path = os.path.join(work_dir, 'chunk_{}_{}.vcf'.format(path_name, chunk_offset))
    job.fileStore.readGlobalFile(vcf_id, vcf_path + '.us')
    
    # Sort the output
    sort_vcf(job, context.runner, vcf_path + '.us', vcf_path)
    
    # do the vcf clip
    left_clip = 0 if chunk_i == 0 else context.config.overlap / 2
    right_clip = 0 if chunk_i == num_chunks - 1 else context.config.overlap / 2
    clip_path = os.path.join(work_dir, 'chunk_{}_{}_clip.vcf'.format(path_name, chunk_offset))
    with open(clip_path, "w") as clip_path_stream:
        offset = vcf_offset + 1
        command=['bcftools', 'view', '-t', '{}:{}-{}'.format(
            path_name, offset + clipped_chunk_offset + left_clip,
            offset + clipped_chunk_offset + context.config.call_chunk_size - right_clip - 1),
                 os.path.basename(vcf_path)]
        context.runner.call(job, command, work_dir=work_dir, outfile=clip_path_stream)

    # save clip.vcf files to job store
    clip_file_id = context.write_intermediate_file(job, clip_path)
    
    return clip_file_id

def run_all_calling(job, context, xg_file_id, chr_gam_ids, chroms, vcf_offsets, sample_name,
                    genotype=False, augment=True, out_name=None):
    """
    Call all the chromosomes and return a merged up vcf/tbi pair
    """
    # we make a child job so that all calling is encapsulated in a top-level job
    child_job = Job()
    job.addChild(child_job)
    vcf_tbi_file_id_pair_list = []
    call_timers_lists = []
    assert len(chr_gam_ids) > 0
    for i in range(len(chr_gam_ids)):
        alignment_file_id = chr_gam_ids[i]
        if len(chr_gam_ids) > 1:
            # 1 gam per chromosome
            chr_label = [chroms[i]]
            chr_offset = [vcf_offsets[i]] if vcf_offsets else [0]
        else:
            # single gam with one or more chromosomes
            chr_label = chroms
            chr_offset = vcf_offsets if vcf_offsets else [0] * len(chroms)
        call_job = child_job.addChildJobFn(run_calling, context, xg_file_id,
                                           alignment_file_id, chr_label, chr_offset,
                                           sample_name, genotype, augment,
                                           cores=context.config.call_chunk_cores,
                                           memory=context.config.call_chunk_mem,
                                           disk=context.config.call_chunk_disk)
        vcf_tbi_file_id_pair_list.append((call_job.rv(0), call_job.rv(1)))
        call_timers_lists.append(call_job.rv(2))
        
    if not out_name:
        out_name = sample_name
    return child_job.addFollowOnJobFn(run_merge_vcf, context, out_name, vcf_tbi_file_id_pair_list,
                                      call_timers_lists,
                                      cores=context.config.call_chunk_cores,
                                      memory=context.config.call_chunk_mem,
                                      disk=context.config.call_chunk_disk).rv()

def run_merge_vcf(job, context, out_name, vcf_tbi_file_id_pair_list, call_timers_lists = []):
    """ Merge up a bunch of chromosome VCFs """

    RealtimeLogger.info("Completed gam merging and gam path variant calling.")
    RealtimeLogger.info("Starting vcf merging vcf files.")

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    timer = TimeTracker('merge-vcf')
    
    vcf_merging_file_key_list = [] 
    for i, vcf_tbi_file_id_pair in enumerate(vcf_tbi_file_id_pair_list):
        vcf_file = os.path.join(work_dir, 'vcf_chunk_{}.vcf.gz'.format(i))
        vcf_file_idx = '{}.tbi'.format(vcf_file)
        job.fileStore.readGlobalFile(vcf_tbi_file_id_pair[0], vcf_file)
        job.fileStore.readGlobalFile(vcf_tbi_file_id_pair[1], vcf_file_idx)
        vcf_merging_file_key_list.append(os.path.basename(vcf_file))

    vcf_merged_file_key = "" 
    if len(vcf_merging_file_key_list) > 1:
        # merge vcf files
        vcf_merged_file_key = "{}.vcf.gz".format(out_name)
        command = ['bcftools', 'concat', '-O', 'z', '-o', os.path.basename(vcf_merged_file_key)]
        command +=  vcf_merging_file_key_list
        context.runner.call(job, command, work_dir=work_dir)
        command=['bcftools', 'tabix', '-f', '-p', 'vcf', os.path.basename(vcf_merged_file_key)]
        context.runner.call(job, command, work_dir=work_dir)
    else:
        vcf_merged_file_key = vcf_merging_file_key_list[0]

    # save variant calling results to the output store
    out_store_key = "{}.vcf.gz".format(out_name)
    vcf_file = os.path.join(work_dir, vcf_merged_file_key)
    vcf_file_idx = vcf_file + ".tbi"
    
    vcf_file_id = context.write_output_file(job, vcf_file, out_store_path = out_store_key)
    vcf_idx_file_id = context.write_output_file(job, vcf_file_idx,
                                                out_store_path = out_store_key + '.tbi')

    # reduce all the timers here from the list of lists
    timer.stop()
    for call_timers in call_timers_lists:
        for call_timer in call_timers:
            timer.add(call_timer)

    if call_timers_lists:
        return vcf_file_id, vcf_idx_file_id, timer
    else:
        return vcf_file_id, vcf_idx_file_id


def run_calling(job, context, xg_file_id, alignment_file_id, path_names, vcf_offsets, sample_name,
                genotype, augment):
    """
    Call a single GAM.  Takes care of splitting the input into chunks based on one or more path,
    processing each chunk in parallel, then merging the result into a single vcf which is returned.
    """
    RealtimeLogger.info("Running variant calling on path(s) {} from alignment file {}".format(','.join(path_names), str(alignment_file_id)))
        
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Tame for work files
    tag = path_names[0] if len(path_names) == 1 else 'chroms'

    # Download the input from the store
    xg_path = os.path.join(work_dir, 'graph.vg.xg')
    job.fileStore.readGlobalFile(xg_file_id, xg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(sample_name))
    job.fileStore.readGlobalFile(alignment_file_id, gam_path)

    # index on coordinates (sort)
    gam_index_path = gam_path + '.index'    
    index_cmd = ['vg', 'index', '-a', os.path.basename(gam_path),
                 '-d', os.path.basename(gam_index_path), '-t', str(context.config.gam_index_cores)]
    timer = TimeTracker('call-gam-index')
    context.runner.call(job, index_cmd, work_dir = work_dir)
    timer.stop()

    # Write a list of paths
    path_list = os.path.join(work_dir, 'path_list.txt')
    offset_map = dict()
    with open(path_list, 'w') as path_list_file:
        for i, path_name in enumerate(path_names):
            path_list_file.write(path_name + '\n')
            offset_map[path_name] = int(vcf_offsets[i]) if vcf_offsets else 0

    # Chunk the graph and gam, using the xg and rocksdb indexes
    output_bed_chunks_path = os.path.join(work_dir, 'output_bed_chunks_{}.bed'.format(tag))
    chunk_cmd = ['vg', 'chunk', '-x', os.path.basename(xg_path),
                 '-a', os.path.basename(gam_index_path), '-c', str(context.config.chunk_context),
                 '-P', os.path.basename(path_list),
                 '-g',
                 '-s', str(context.config.call_chunk_size),
                 '-o', str(context.config.overlap),
                 '-b', 'call_chunk_{}'.format(tag),
                 '-t', str(context.config.call_chunk_cores),
                 '-E', os.path.basename(output_bed_chunks_path),
                 '-f']
    timer.start('call-chunk')
    context.runner.call(job, chunk_cmd, work_dir=work_dir)
    timer.stop()

    # Scrape the BED into memory
    bed_lines = []
    path_bounds = dict()    
    with open(output_bed_chunks_path) as output_bed:
        for line in output_bed:
            toks = line.split('\t')
            if len(toks) > 3:
                bed_lines.append(toks)
                chrom, start, end = toks[0], int(toks[1]), int(toks[2])
                if chrom not in path_bounds:
                    path_bounds[chrom] = (start, end)
                else:
                    path_bounds[chrom] = (min(start, path_bounds[chrom][0]),
                                              max(end, path_bounds[chrom][1]))

    # Infer the size of the path from our BED (a bit hacky)
    path_size = dict()
    for name, bounds in path_bounds.items():
        path_size[name] = path_bounds[name][1] - path_bounds[name][0]

    # Keep track of offset in each path
    cur_path_offset = defaultdict(int)

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)
        
    # Go through the BED output of vg chunk, adding a child calling job for
    # each chunk                
    clip_file_ids = []
    call_timers = [timer]
    for toks in bed_lines:
        chunk_bed_chrom = toks[0]
        chunk_bed_start = int(toks[1])
        chunk_bed_end = int(toks[2])
        chunk_bed_size = chunk_bed_end - chunk_bed_start
        gam_chunk_path = os.path.join(work_dir, os.path.basename(toks[3].strip()))
        assert os.path.splitext(gam_chunk_path)[1] == '.gam'
        vg_chunk_path = os.path.splitext(gam_chunk_path)[0] + '.vg'
        gam_chunk_file_id = context.write_intermediate_file(job, gam_chunk_path)
        vg_chunk_file_id = context.write_intermediate_file(job, vg_chunk_path)
        chunk_i = cur_path_offset[chunk_bed_chrom]        
        clipped_chunk_offset = chunk_i * context.config.call_chunk_size - chunk_i * context.config.overlap
        cur_path_offset[chunk_bed_chrom] += 1

        call_job =  child_job.addChildJobFn(run_call_chunk, context, chunk_bed_chrom, chunk_i,
                                            len(bed_lines),
                                            chunk_bed_start, clipped_chunk_offset,
                                            None, vg_chunk_file_id, gam_chunk_file_id,
                                            path_size[chunk_bed_chrom], offset_map[chunk_bed_chrom],
                                            sample_name, genotype, augment,
                                            cores=context.config.misc_cores,
                                            memory=context.config.misc_mem, disk=context.config.misc_disk)
        clip_file_ids.append(call_job.rv(0))
        call_timers.append(call_job.rv(1))
        
    merge_job = child_job.addFollowOnJobFn(run_merge_vcf_chunks, context, tag,
                                           clip_file_ids,
                                           cores=context.config.call_chunk_cores,
                                           memory=context.config.call_chunk_mem,
                                           disk=context.config.call_chunk_disk)
        
    vcf_out_file_id = merge_job.rv(0)
    tbi_out_file_id = merge_job.rv(1)
    
    return vcf_out_file_id, tbi_out_file_id, call_timers


def run_merge_vcf_chunks(job, context, path_name, clip_file_ids):
    """ merge a bunch of clipped vcfs created above, taking care to 
    fix up the headers.  everything expected to be sorted already """
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vcf_path = os.path.join(work_dir, path_name + ".vcf")
    
    for chunk_i, clip_file_id in enumerate(clip_file_ids):
        
        # Download clip.vcf file from the store
        clip_path = os.path.join(work_dir, 'clip_{}.vcf'.format(chunk_i))
        job.fileStore.readGlobalFile(clip_file_id, clip_path)

        if chunk_i == 0:
            # copy everything including the header
            with open(vcf_path, "w") as outfile:
                context.runner.call(job, ['cat', os.path.basename(clip_path)], outfile=outfile,
                                     work_dir=work_dir)
        else:
            # add on everythin but header
            with open(vcf_path, "a") as outfile:
                context.runner.call(job, ['bcftools', 'view', '-H', os.path.basename(clip_path)],
                                     outfile=outfile, work_dir=work_dir)

    # add a compressed indexed version
    vcf_gz_file = vcf_path + ".gz"
    with open(vcf_gz_file, "w") as vcf_gz_file_stream:
        command=['bgzip', '-c', '{}'.format(os.path.basename(vcf_path))]
        context.runner.call(job, command, work_dir=work_dir, outfile=vcf_gz_file_stream)
    command=['bcftools', 'tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
    context.runner.call(job, command, work_dir=work_dir)

    # Save merged vcf files to the job store
    vcf_gz_file_id = context.write_intermediate_file(job, vcf_path+".gz")
    vcf_tbi_file_id = context.write_intermediate_file(job, vcf_path+".gz.tbi")

    RealtimeLogger.info("Completed variant calling on path {}".format(path_name))

    return vcf_gz_file_id, vcf_tbi_file_id
    
def call_main(context, options):
    """ entrypoint for calling """

    validate_call_options(options)
            
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()

            # Upload local files to the job store
            inputXGFileID = toil.importFile(options.xg_path)
            inputGamFileIDs = []
            for inputGamFileID in options.gams:
                inputGamFileIDs.append(toil.importFile(inputGamFileID))

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_all_calling, context, inputXGFileID, inputGamFileIDs,
                                     options.chroms, options.vcf_offsets, options.sample_name,
                                     options.genotype, not options.no_augment,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)
            init_job.addFollowOn(root_job)            
            
            # Run the job and store the returned list of output files to download
            toil.start(init_job)
        else:
            toil.restart()
                
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
    
