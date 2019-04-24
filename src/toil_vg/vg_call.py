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
                        " Must be same length/order as --gams.  All paths used if not specified")
    # todo: move to chunked_call_parse_args and share with toil-vg run
    parser.add_argument("--gams", nargs='+', required=True, type=make_url,
                        help="GAMs to call.  One per chromosome in the same order as --chroms, or just one. "
                        " Indexes (.gai) will be used if found.")
    
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
                        help="chunk size (set to 0 to disable chunking)")
    parser.add_argument("--call_opts", type=str,
                        help="arguments to pass to vg call (wrapped in \"\")")
    parser.add_argument("--genotype", action="store_true",
                        help="use vg genotype instead of vg call")
    parser.add_argument("--genotype_opts", type=str,
                        help="arguments to pass to vg genotype (wrapped in \"\")")
    parser.add_argument("--filter_opts", type=str,
                        help="argument to pass to vg filter (wrapped in \"\")")
    parser.add_argument("--calling_cores", type=int,
                        help="number of threads during the variant calling step")
    parser.add_argument("--calling_mem", type=str,
                        help="memory alotment during the variant calling step")
    parser.add_argument("--call_chunk_cores", type=int,
                        help="number of threads used for extracting chunks for calling")
    parser.add_argument("--call_chunk_mem", type=str,
                        help="memory alotment for extracting chunks for calling")
    parser.add_argument("--vcf_offsets", nargs='+', default=[],
                        help="offset(s) to apply to output vcfs(s). (order of --chroms)")
    parser.add_argument("--recall", action="store_true",
                        help="only call variants present in the graph")
    parser.add_argument("--alt_path_gam", type=make_url,
                        help="GAM storing alt paths for use with --genotype_vcf, as created by toil-vg construct/index"
                        " --alt_path_gam_index. Must have associated .gai")
    parser.add_argument("--genotype_vcf", type=make_url,
                        help="Genotype variants in a given .vcf.gz file. This file must have been used to construct"
                        " the graph (with --alt_path_gam_index")

def validate_call_options(options):    
    require(not options.chroms or len(options.chroms) == len(options.gams) or len(options.gams) == 1,
            'Number of --chroms must be 1 or same as number of --gams')
    require(not options.vcf_offsets or len(options.vcf_offsets) == len(options.chroms),
            'Number of --vcf_offsets if specified must be same as number of --chroms')
    require(not options.genotype or not options.recall,
            '--recall not yet supported with --genotype')
    require((options.alt_path_gam is not None) == (options.genotype_vcf is not None),
            '--alt_path_gam must be used in conjunction with --genotype_vcf')
    require(not options.genotype_vcf or options.genotype_vcf.endswith('.vcf.gz'),
            'file passed to --genotype_vcf must end with .vcf.gz')
    
def sort_vcf(job, drunner, vcf_path, sorted_vcf_path):
    """ from vcflib """
    vcf_dir, vcf_name = os.path.split(vcf_path)
    with open(sorted_vcf_path, "w") as outfile:
        drunner.call(job, [['bcftools', 'view', '-h', vcf_name]], outfile=outfile,
                     work_dir=vcf_dir)
        drunner.call(job, [['bcftools', 'view', '-H', vcf_name],
                      ['sort', '-k1,1d', '-k2,2n']], outfile=outfile,
                     work_dir=vcf_dir)
        
def run_vg_call(job, context, sample_name, vg_id, gam_id, xg_id = None,
                path_names = [], seq_names = [], seq_offsets = [], seq_lengths = [],
                chunk_name = 'call', genotype = False, recall = False, clip_info = None,
                augment_results = None, augment_only = False, alt_gam_id = None,
                genotype_vcf_id = None, genotype_tbi_id = None):
    """ Run vg call or vg genotype on a single graph.

    Returns (vcf_id, pileup_id, xg_id, gam_id, augmented_graph_id). pileup_id and xg_id
    can be same as input if they are not computed.  If pileup/xg/augmented are 
    computed, the returned ids will be None unless appropriate keep_flag set
    (to prevent sending them to the file store if they aren't wanted)

    User is responsible to make sure that options passed in context.config.*_opts don't conflict
    with seq_names, seq_offsets, seq_lengths etc. If not provided, the pileup is computed.

    gam filtering is only done if filter_opts are passed in. 

    chunk_name option is only for working filenames (to make more readable)
    
    When running vg genotype, we can't recall (since recall in vg genotype needs a VCF). 

    augment_results is a dict with the ids of augment results from a previous run of augment
    if it's given, then they are used and augment is not run

    if augment_only is True, then calling is skipped and augment_results will be output

    """
    
    work_dir = job.fileStore.getLocalTempDir()

    filter_opts = context.config.filter_opts if not recall else context.config.recall_filter_opts
    augment_opts = context.config.augment_opts
    if genotype:
        call_opts = context.config.genotype_opts
    else:
        call_opts = context.config.call_opts if not recall else context.config.recall_opts

    # Read our input files from the store
    vg_path = os.path.join(work_dir, '{}.vg'.format(chunk_name))
    job.fileStore.readGlobalFile(vg_id, vg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(chunk_name))
    if not augment_results:
        job.fileStore.readGlobalFile(gam_id, gam_path)
        xg_path = os.path.join(work_dir, '{}.xg'.format(chunk_name))
        defray = filter_opts and ('-D' in filter_opts or '--defray-ends' in filter_opts)
        if xg_id and defray:
            job.fileStore.readGlobalFile(xg_id, xg_path)
        if alt_gam_id:
            alt_gam_path = os.path.join(work_dir, '{}_alts.gam'.format(chunk_name))
            job.fileStore.readGlobalFile(alt_gam_id, alt_gam_path)
                                        
        
    # Define paths for all the files we might make
    pu_path = os.path.join(work_dir, '{}.pu'.format(chunk_name))
    trans_path = os.path.join(work_dir, '{}.trans'.format(chunk_name))
    support_path = os.path.join(work_dir, '{}.support'.format(chunk_name))
    aug_path = os.path.join(work_dir, '{}_aug.vg'.format(chunk_name))
    aug_gam_path = os.path.join(work_dir, '{}_aug.gam'.format(chunk_name))
    vcf_path = os.path.join(work_dir, '{}_call.vcf'.format(chunk_name))
    sorted_vcf_path = os.path.join(work_dir, '{}_call_sorted.vcf'.format(chunk_name))

    timer = TimeTracker()

    if not augment_results:
        # we only need an xg if using vg filter -D
        if not xg_id and defray:
            timer.start('chunk-xg')
            context.runner.call(job, ['vg', 'index', os.path.basename(vg_path), '-x',
                                       os.path.basename(xg_path), '-t', str(context.config.calling_cores)],
                                 work_dir = work_dir)
            timer.stop()

        # optional alt path augmentation
        if alt_gam_id:
            # Note: There's no reason why this should modify the nodes and edges as the paths were already
            # embedded.  But if it does due to a bug somewhere, terrible downstream errors will result. 
            alt_augment_cmd = ['vg', 'augment', '-i', os.path.basename(vg_path), os.path.basename(alt_gam_path)]
            vg_alt_aug_path = os.path.join(work_dir, '{}_alts.vg'.format(chunk_name))
            with open(vg_alt_aug_path, 'w') as aug_vg_file:
                context.runner.call(job, alt_augment_cmd, work_dir=work_dir, outfile=aug_vg_file)
            vg_path = vg_alt_aug_path

        # optional gam filtering
        gam_filter_path = gam_path + '.filter'
        filter_command = None
        if filter_opts:
            filter_command = ['vg', 'filter', os.path.basename(gam_path), '-t', '1'] + filter_opts
            if defray:
                filter_command += ['-x', os.path.basename(xg_path)]
            if genotype:
                with open(gam_filter_path, 'w') as gam_filter_file:
                    context.runner.call(job, filter_command, work_dir = work_dir, outfile = gam_filter_file)
                filter_command = None
        else:
            gam_filter_path = gam_path

        # augmentation
        augment_command = []
        augment_generated_opts = []
        if filter_command is not None:
            aug_gam_input = '-'
            augment_command.append(filter_command)
        else:
            aug_gam_input = os.path.basename(gam_filter_path)

        if genotype:
            augment_generated_opts = ['--augmentation-mode', 'direct',
                                      '--alignment-out', os.path.basename(aug_gam_path)]
            augment_opts = []
        else:
            augment_generated_opts = ['--augmentation-mode', 'pileup',
                                      '--translation', os.path.basename(trans_path),
                                      '--support', os.path.basename(support_path)]
            if recall:
                augment_opts = []
                augment_generated_opts += ['--recall']

        augment_command.append(['vg', 'augment', os.path.basename(vg_path), aug_gam_input,
                        '-t', str(context.config.calling_cores)] + augment_opts + augment_generated_opts)

        try:
            with open(aug_path, 'w') as aug_stream:
                timer.start('call-filter-augment')
                context.runner.call(job, augment_command, work_dir=work_dir, outfile=aug_stream)
                timer.stop()
        except Exception as e:
            logging.error("Augmentation failed. Dumping files.")
            for dump_path in [vg_path, gam_path, gam_filter_path]:
                if dump_path and os.path.isfile(dump_path):
                    context.write_output_file(job, dump_path)        
            raise e
    else:
        # we download the augment output instead of running it
        job.fileStore.readGlobalFile(augment_results['aug-graph'], aug_path)
        if not genotype:
            job.fileStore.readGlobalFile(augment_results['support'], support_path)
            job.fileStore.readGlobalFile(augment_results['translation'], trans_path)

    # We're going to stop here and return our augmentation results
    if augment_only:
        augment_output_results = { 'aug-graph' : context.write_intermediate_file(job, aug_path) }
        if not genotype:
            augment_output_results['support'] = context.write_intermediate_file(job, support_path)
            augment_output_results['translation'] = context.write_intermediate_file(job, trans_path)
        return augment_output_results

    # naming options shared between call and genotype
    name_opts = []
    for path_name in path_names:
        name_opts += ['-r', path_name]
    for seq_name in seq_names:
        name_opts += ['-c', seq_name]
    for seq_length in seq_lengths:
        name_opts += ['-l', seq_length]
    if not genotype_vcf_id:
        for seq_offset in seq_offsets:
            name_opts += ['-o', seq_offset]
                    
    if genotype:
        # We can't do recall (no passed VCF) 
        assert(not recall)
        
        # How do we actually genotype
        # (Genotype runs fine with its built-in augmentation, but I suspect it causes trouble
        # with some CI tests, perhaps by taking more memory, so we leave augmentation in its
        # own command)
        command = ['vg', 'genotype', os.path.basename(aug_path), '-t',
                   str(context.config.calling_cores), '-s', sample_name,
                   '-v', '-E', '-G', os.path.basename(aug_gam_path)]
                           
        if call_opts:
            command += call_opts
        command += name_opts
        
        try:
            with open(vcf_path, 'w') as vggenotype_stdout:
                timer.start('genotype')
                context.runner.call(job, command, work_dir=work_dir,
                                     outfile=vggenotype_stdout)
                timer.stop()

        except Exception as e:
            logging.error("Genotyping failed. Dumping files.")
            for dump_path in [vg_path, aug_gam_path, aug_path]:
                if dump_path and os.path.isfile(dump_path):
                    context.write_output_file(job, dump_path)        
            raise e
                
    else:
        # call
        command = ['vg', 'call', os.path.basename(aug_path), '-t',
                   str(context.config.calling_cores), '-S', sample_name,
                   '-z', os.path.basename(trans_path),
                   '-s', os.path.basename(support_path),
                   '-b', os.path.basename(vg_path)]
                
        if call_opts:
            command += call_opts
        command += name_opts

        if genotype_vcf_id:
            genotype_vcf_path = os.path.join(work_dir, '{}_to_genotype.vcf.gz'.format(chunk_name))
            job.fileStore.readGlobalFile(genotype_vcf_id, genotype_vcf_path)
            job.fileStore.readGlobalFile(genotype_tbi_id, genotype_vcf_path + '.tbi')
            if clip_info and context.config.call_chunk_size != 0:
                genotype_clip_vcf_path = os.path.join(work_dir, '{}_to_genotype_clipped.vcf.gz'.format(chunk_name))
                # todo: share code with clipping below
                left_clip = 0 if clip_info['chunk_i'] == 0 else context.config.overlap / 2
                right_clip = 0 if clip_info['chunk_i'] == clip_info['chunk_n'] - 1 else context.config.overlap / 2
                offset = clip_info['offset'] + 1 # adjust to 1-based vcf
                clip_command=['bcftools', 'view', '-O', 'z', '-t', '{}:{}-{}'.format(
                    path_name, offset + clip_info['clipped_chunk_offset'] + left_clip,
                    offset + clip_info['clipped_chunk_offset'] + context.config.call_chunk_size - right_clip - 1),
                         os.path.basename(genotype_vcf_path), '--output-file', os.path.basename(genotype_clip_vcf_path)]
                context.runner.call(job, clip_command, work_dir=work_dir)
                context.runner.call(job, ['tabix', '-f', '-p', 'vcf', os.path.basename(genotype_clip_vcf_path)], work_dir=work_dir)
                genotype_vcf_path = genotype_clip_vcf_path
            command += ['-f', os.path.basename(genotype_vcf_path)]

        try:
            with open(vcf_path, 'w') as vgcall_stdout:
                timer.start('call')
                context.runner.call(job, command, work_dir=work_dir, outfile=vgcall_stdout)
                timer.stop()            
        except Exception as e:
            logging.error("Calling failed. Dumping files.")
            for dump_path in [vg_path, pu_path,
                              aug_path, support_path, trans_path, aug_gam_path]:
                if dump_path and os.path.isfile(dump_path):
                    context.write_output_file(job, dump_path)
            if genotype_vcf_id and os.path.isfile(genotype_vcf_path):
                context.write_output_file(job, genotype_vcf_path)
                context.write_output_file(job, genotype_vcf_path + '.tbi')
            raise e

    # Sort the output
    sort_vcf(job, context.runner, vcf_path, sorted_vcf_path)

    # Optional clip
    if clip_info and context.config.call_chunk_size != 0 and not genotype_vcf_id:
        left_clip = 0 if clip_info['chunk_i'] == 0 else context.config.overlap / 2
        right_clip = 0 if clip_info['chunk_i'] == clip_info['chunk_n'] - 1 else context.config.overlap / 2
        clip_path = os.path.join(work_dir, '{}_call_clip.vcf'.format(chunk_name))
        offset = clip_info['offset'] + 1 # adjust to 1-based vcf
        command=['bcftools', 'view', '-t', '{}:{}-{}'.format(
            path_name, offset + clip_info['clipped_chunk_offset'] + left_clip,
            offset + clip_info['clipped_chunk_offset'] + context.config.call_chunk_size - right_clip - 1),
                 os.path.basename(sorted_vcf_path), '--output-file', os.path.basename(clip_path)]
        context.runner.call(job, command, work_dir=work_dir)
        vcf_id = context.write_intermediate_file(job, clip_path)
    else:
        vcf_id =  context.write_intermediate_file(job, sorted_vcf_path)
        
    return vcf_id, timer

# Get the lengths of our paths in the graph
def run_xg_paths(job, context, xg_id):
    """ pull path names and sizes from the xg """
    work_dir = job.fileStore.getLocalTempDir()
    xg_path = os.path.join(work_dir, 'index.xg')
    job.fileStore.readGlobalFile(xg_id, xg_path)    
    vg_paths_output = context.runner.call(job, ['vg', 'paths', '--xg', os.path.basename(xg_path), '--lengths'],
                                          work_dir = work_dir, check_output = True)
    path_size = dict([line.strip().split('\t') for line in vg_paths_output.split('\n') if len(line) > 2])
    return path_size

def run_all_calling(job, context, xg_file_id, chr_gam_ids, chr_gam_idx_ids, chroms, vcf_offsets, sample_name,
                    genotype=False, out_name=None, recall=False, alt_gam_id=None, alt_gai_id=None,
                    genotype_vcf_id=None, genotype_tbi_id=None):
    path_sizes_job = job.addChildJobFn(run_xg_paths, context, xg_file_id,
                                       memory=context.config.call_chunk_mem,
                                       disk=context.config.call_chunk_disk)
    path_sizes = path_sizes_job.rv()
    calling_job = path_sizes_job.addFollowOnJobFn(run_all_calling2, context, xg_file_id, chr_gam_ids, chr_gam_idx_ids,
                                                  chroms, path_sizes,  vcf_offsets, sample_name, genotype, out_name,
                                                  recall, alt_gam_id, alt_gai_id, genotype_vcf_id, genotype_tbi_id)
    return calling_job.rv()

def run_all_calling2(job, context, xg_file_id, chr_gam_ids, chr_gam_idx_ids, chroms, path_sizes, vcf_offsets, sample_name,
                     genotype=False, out_name=None, recall=False, alt_gam_id=None, alt_gai_id=None,
                     genotype_vcf_id=None, genotype_tbi_id=None):
    """
    Call all the chromosomes and return a merged up vcf/tbi pair
    """
    # we make a child job so that all calling is encapsulated in a top-level job
    child_job = Job()
    job.addChild(child_job)
    vcf_ids = []
    tbi_ids = []
    call_timers_lists = []
    assert len(chr_gam_ids) > 0
    if not chr_gam_idx_ids:
        chr_gam_idx_ids = [None] * len(chr_gam_ids)
    if not chroms:
        chroms = [name for name in path_sizes.keys() if path_sizes[name] > 0]
    assert len(chr_gam_ids) == len(chr_gam_idx_ids)
    for i in range(len(chr_gam_ids)):
        alignment_file_id = chr_gam_ids[i]
        alignment_index_id = chr_gam_idx_ids[i]
        if len(chr_gam_ids) > 1:
            # 1 gam per chromosome
            chr_label = [chroms[i]]
            chr_offset = [vcf_offsets[i]] if vcf_offsets else [0]
        else:
            # single gam with one or more chromosomes
            chr_label = chroms
            chr_offset = vcf_offsets if vcf_offsets else [0] * len(chroms)
        chunk_job = child_job.addChildJobFn(run_chunking, context, xg_file_id,
                                            alignment_file_id, alignment_index_id, chr_label, chr_offset, path_sizes,
                                            sample_name, genotype=genotype, recall=recall,
                                            alt_gam_id=alt_gam_id, alt_gai_id=alt_gai_id,
                                            genotype_vcf_id=genotype_vcf_id,
                                            genotype_tbi_id=genotype_tbi_id,
                                            cores=context.config.call_chunk_cores,
                                            memory=context.config.call_chunk_mem,
                                            disk=context.config.call_chunk_disk)
        call_job = chunk_job.addFollowOnJobFn(run_chunked_calling, context, chunk_job.rv(0),
                                              genotype, recall, chunk_job.rv(1),
                                              cores=context.config.misc_cores,
                                              memory=context.config.misc_mem,
                                              disk=context.config.misc_disk)
        vcf_ids.append(call_job.rv(0))
        tbi_ids.append(call_job.rv(1))
        call_timers_lists.append(call_job.rv(2))
        
    if not out_name:
        out_name = sample_name
    return child_job.addFollowOnJobFn(run_concat_vcfs, context, out_name, vcf_ids, tbi_ids,
                                      write_to_outstore = True,
                                      call_timers_lists = call_timers_lists,
                                      cores=context.config.call_chunk_cores,
                                      memory=context.config.call_chunk_mem,
                                      disk=context.config.call_chunk_disk).rv()

def run_concat_vcfs(job, context, out_name, vcf_ids, tbi_ids = [], write_to_outstore = False,
                    call_timers_lists = []):
    """ Concat up a bunch of VCFs. Input assumed to be bgzipped iff tbi_ids specified """

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    timer = TimeTracker('merge-vcf')

    # Download the input
    vcf_paths = []
    tbi_paths = []
    for i, vcf_id in enumerate(vcf_ids):
        vcf_path = os.path.join(work_dir, 'vcf_chunk_{}.vcf'.format(i))
        if len(tbi_ids) == len(vcf_ids):
            vcf_path += '.gz'
            tbi_path = vcf_path + '.tbi'
            job.fileStore.readGlobalFile(tbi_ids[i], tbi_path)
            tbi_paths.append(tbi_path)
        job.fileStore.readGlobalFile(vcf_id, vcf_path)
        vcf_paths.append(vcf_path)

    out_file = os.path.join(work_dir, out_name + '.vcf.gz')

    # Merge with bcftools
    command = ['bcftools', 'concat'] + [os.path.basename(vcf_path) for vcf_path in vcf_paths]
    command += ['--output', os.path.basename(out_file), '--output-type', 'z']
    context.runner.call(job, command, work_dir=work_dir)
    context.runner.call(job, ['tabix', '--force', '--preset', 'vcf', os.path.basename(out_file)],
                        work_dir = work_dir)

    if write_to_outstore:
        vcf_file_id = context.write_output_file(job, out_file)
        vcf_idx_file_id = context.write_output_file(job, out_file +'.tbi')
    else:
        vcf_file_id = context.write_intermediate_file(job, out_file)
        vcf_idx_file_id = context.write_intermediate_file(job, out_file + '.tbi')

    # reduce all the timers here from the list of lists
    timer.stop()
    for call_timers in call_timers_lists:
        for call_timer in call_timers:
            timer.add(call_timer)

    if call_timers_lists:
        return vcf_file_id, vcf_idx_file_id, timer
    else:
        return vcf_file_id, vcf_idx_file_id

def run_chunking(job, context, xg_file_id, alignment_file_id, alignment_index_id, path_names, vcf_offsets, path_sizes, sample_name,
                 genotype, recall, alt_gam_id, alt_gai_id, genotype_vcf_id, genotype_tbi_id):
    """
    Split a gam and xg up into a bunch of vg/gam pairs, one for each calling chunk.  also keep track of the 
    various offsets needed to call and clip them.
    """

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Name for work files
    tag = path_names[0] if len(path_names) == 1 else 'chroms'

    # Download the input from the store
    xg_path = os.path.join(work_dir, 'graph.xg')
    job.fileStore.readGlobalFile(xg_file_id, xg_path)
    gam_path = os.path.join(work_dir, '{}_{}.gam'.format(sample_name, tag))
    job.fileStore.readGlobalFile(alignment_file_id, gam_path)
    alt_gam_path = os.path.join(work_dir, '{}_{}_alts.gam'.format(sample_name, tag))
    if alt_gam_id:
        job.fileStore.readGlobalFile(alt_gam_id, alt_gam_path)
        job.fileStore.readGlobalFile(alt_gai_id, alt_gam_path + '.gai')
    
    # This will be a list of dicts, each one corresponding to the info for a chunk
    output_chunk_info = []
    
    # Bypass chunking when chunk_size set to 0
    if context.config.call_chunk_size == 0:
        timer = TimeTracker('call-chunk-bypass')
        # convert the xg to vg
        context.runner.call(job, ['vg', 'xg', '-i', os.path.basename(xg_path), '-X', 'graph.vg'],
                            work_dir = work_dir)
        vg_id = context.write_intermediate_file(job, os.path.join(work_dir, 'graph.vg'))
        timer.stop()

        # return a job for each path so we can run them in parallel, but they will all
        # use the same graph and gam. 
        for chunk_i, path_name in enumerate(path_names):
            chunk_info = {
                'chrom' : path_name,
                'chunk_i' : 0,
                'chunk_n' : 1,
                'chunk_start' : 0,
                'clipped_chunk_offset' : 0,
                'vg_id' : vg_id,
                'gam_id' : alignment_file_id,
                'xg_id' : xg_file_id,
                'path_size' : path_sizes[path_name],
                'offset' : 0 if not vcf_offsets or chunk_i >= len(path_names) else vcf_offsets[chunk_i],
                'sample' : sample_name,
                'alt_gam_id' : alt_gam_id,
                'genotype_vcf_id' : genotype_vcf_id,
                'genotype_tbi_id' : genotype_tbi_id}
            output_chunk_info.append(chunk_info)
            
        return output_chunk_info, [timer]

    # Sort and index the GAM file if index not provided
    timer = TimeTracker('call-gam-index')    
    if alignment_index_id:
        gam_sort_path = gam_path
        gam_index_path = gam_sort_path + '.gai'
        job.fileStore.readGlobalFile(alignment_index_id, gam_index_path)
    else:
        gam_sort_path = gam_path + '.sorted.gam'
        gam_index_path = gam_sort_path + '.gai'
        with open(gam_sort_path, "w") as gam_sort_stream:
            gamsort_cmd = ['vg', 'gamsort', '-i', os.path.basename(gam_index_path), os.path.basename(gam_path),
                           '--threads', str(context.config.call_chunk_cores)]
            context.runner.call(job, gamsort_cmd, work_dir=work_dir, outfile=gam_sort_stream)
        # We may spend a lot of time on these, so drop them in the output store.
        if context.config.keep_sorted_gams:
            context.write_output_file(job, gam_sort_path)
            context.write_output_file(job, gam_index_path)
    timer.stop()
    
    # Write a list of paths
    path_list = os.path.join(work_dir, 'path_list.txt')
    offset_map = dict()
    with open(path_list, 'w') as path_list_file:
        for i, path_name in enumerate(path_names):
            path_list_file.write(path_name + '\n')
            offset_map[path_name] = int(vcf_offsets[i]) if vcf_offsets else 0

    # Apply chunk override for recall.
    # Todo: fix vg chunk to only expand variants (not reference) so that we can leave this super high
    # all the time
    
    context_size = int(context.config.chunk_context) if not recall else int(context.config.recall_context)
    
    # Chunk the graph and gam, using the xg and rocksdb indexes.
    # GAM index isn't passed but it needs to be next to the GAM file.
    output_bed_chunks_path = os.path.join(work_dir, 'output_bed_chunks_{}.bed'.format(tag))
    chunk_prefix = 'call_chunk_{}'.format(tag)
    chunk_cmd = ['vg', 'chunk', '-x', os.path.basename(xg_path),
                 '-a', os.path.basename(gam_sort_path), '-c', str(context_size),
                 '-P', os.path.basename(path_list),
                 '-g',
                 '-s', str(context.config.call_chunk_size),
                 '-o', str(context.config.overlap),
                 '-b', chunk_prefix,
                 '-t', str(context.config.call_chunk_cores),
                 '-E', os.path.basename(output_bed_chunks_path),
                 '-f']
    if alt_gam_id:
        chunk_cmd += ['-a', os.path.basename(alt_gam_path)]
    timer.start('call-chunk')
    context.runner.call(job, chunk_cmd, work_dir=work_dir)
    timer.stop()

    # Scrape the BED into memory
    bed_lines = []
    chunk_counts = defaultdict(int)
    with open(output_bed_chunks_path) as output_bed:
        for line in output_bed:
            toks = line.split('\t')
            if len(toks) > 3:
                bed_lines.append(toks)
                chunk_counts[toks[0]] += 1

    # Keep track of offset in each path
    cur_path_offset = defaultdict(int)
        
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
        if alt_gam_id:
            alt_gam_chunk_path = gam_chunk_path.replace(chunk_prefix, chunk_prefix + '-1')
            alt_gam_chunk_file_id = context.write_intermediate_file(job, alt_gam_chunk_path)

        chunk_info = {
            'chrom' : chunk_bed_chrom,
            'chunk_i' : chunk_i,
            'chunk_n' : chunk_counts[chunk_bed_chrom],
            'chunk_start' : chunk_bed_start,
            'clipped_chunk_offset' : clipped_chunk_offset,
            'vg_id' : vg_chunk_file_id,
            'gam_id' : gam_chunk_file_id,
            'xg_id' : None,
            'path_size' : path_sizes[chunk_bed_chrom],
            'offset' : offset_map[chunk_bed_chrom],
            'sample' : sample_name,
            'alt_gam_id' : alt_gam_chunk_file_id if alt_gam_id else None,
            'genotype_vcf_id' : genotype_vcf_id,
            'genotype_tbi_id' : genotype_tbi_id
        }
        output_chunk_info.append(chunk_info)

    return output_chunk_info, call_timers

def run_chunked_calling(job, context, chunk_infos, genotype, recall, call_timers):
    """
    spawn a calling job for each chunk then merge them together
    """
    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)

    path_names = set()

    # If no chunking and many paths, we augment once first and not before calling
    # so we don't waste resources augmenting the same graph again and again
    # Note: should only do this when len(chunk_infos) > 1, but leaving as is so the tests hit it!
    if context.config.call_chunk_size == 0:
        chunk_info = chunk_infos[0]
        augment_job = child_job.addChildJobFn(
            run_vg_call,
            context,
            chunk_info['sample'],
            chunk_info['vg_id'],
            chunk_info['gam_id'],
            xg_id = chunk_info['xg_id'],
            path_names = [chunk_info['chrom']],
            seq_names = [chunk_info['chrom']],
            seq_offsets = [chunk_info['chunk_start'] + chunk_info['offset']],
            seq_lengths = [chunk_info['path_size']],
            chunk_name = 'chunk_{}_{}'.format(chunk_info['chrom'], chunk_info['chunk_start']),
            genotype = genotype,
            recall = recall,
            clip_info = chunk_info,
            augment_only = True,
            alt_gam_id = chunk_info['alt_gam_id'],
            cores=context.config.calling_cores,
            memory=context.config.calling_mem, disk=context.config.calling_disk)
        augment_results = augment_job.rv()
        next_job = Job()
        augment_job.addFollowOn(next_job)
        child_job = next_job
    else:
        augment_results = None
    
    clip_file_ids = []
    for chunk_info in chunk_infos:
        path_names.add(chunk_info['chrom'])

        # Run vg call
        call_job = child_job.addChildJobFn(
            run_vg_call,
            context,
            chunk_info['sample'],
            chunk_info['vg_id'],
            chunk_info['gam_id'],
            xg_id = chunk_info['xg_id'],
            path_names = [chunk_info['chrom']],
            seq_names = [chunk_info['chrom']],
            seq_offsets = [chunk_info['chunk_start'] + chunk_info['offset']],
            seq_lengths = [chunk_info['path_size']],
            chunk_name = 'chunk_{}_{}'.format(chunk_info['chrom'], chunk_info['chunk_start']),
            genotype = genotype,
            recall = recall,
            clip_info = chunk_info,
            alt_gam_id = chunk_info['alt_gam_id'],
            genotype_vcf_id = chunk_info['genotype_vcf_id'],
            genotype_tbi_id = chunk_info['genotype_tbi_id'],
            augment_results = augment_results,
            cores=context.config.calling_cores,
            memory=context.config.calling_mem, disk=context.config.calling_disk)
        vcf_id, call_timer = call_job.rv(0), call_job.rv(1)
        
        clip_file_ids.append(vcf_id)
        call_timers.append(call_timer)

    tag = list(path_names)[0] if len(path_names) == 1 else 'chroms'
        
    merge_job = child_job.addFollowOnJobFn(run_concat_vcfs, context, tag,
                                           clip_file_ids,
                                           cores=context.config.call_chunk_cores,
                                           memory=context.config.call_chunk_mem,
                                           disk=context.config.call_chunk_disk)
        
    vcf_out_file_id = merge_job.rv(0)
    tbi_out_file_id = merge_job.rv(1)
    
    return vcf_out_file_id, tbi_out_file_id, call_timers
    
def call_main(context, options):
    """ entrypoint for calling """

    validate_call_options(options)
            
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)

            # Upload local files to the job store
            inputXGFileID = importer.load(options.xg_path)
            inputGamFileIDs = []
            inputGamIndexFileIDs = []
            for inputGam in options.gams:
                inputGamFileIDs.append(importer.load(inputGam))
                try:
                    inputGamIndexID = toil.importFile(inputGam + '.gai')
                except:
                    inputGamIndexID = None
                # we allow some GAMs to have indexes and some to have None
                inputGamIndexFileIDs.append(inputGamIndexID)
            if options.alt_path_gam is not None:
                inputAltPathGamID = importer.load(options.alt_path_gam)
                inputAltPathGaiID = toil.importFile(options.alt_path_gam + '.gai')
            else:
                inputAltPathGamID, inputAltPathGaiID = None, None
            if options.genotype_vcf is not None:
                inputVcfID = importer.load(options.genotype_vcf)
                inputTbiID = importer.load(options.genotype_vcf + '.tbi')
            else:
                inputVcfID, inputTbiID = None, None

            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_all_calling, context, importer.resolve(inputXGFileID),
                                     importer.resolve(inputGamFileIDs), inputGamIndexFileIDs,
                                     options.chroms, options.vcf_offsets, options.sample_name,
                                     genotype=options.genotype,
                                     recall=options.recall or options.genotype_vcf is not None,
                                     alt_gam_id=importer.resolve(inputAltPathGamID),
                                     alt_gai_id=inputAltPathGaiID,
                                     genotype_vcf_id=importer.resolve(inputVcfID),
                                     genotype_tbi_id=importer.resolve(inputTbiID),
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
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
    
