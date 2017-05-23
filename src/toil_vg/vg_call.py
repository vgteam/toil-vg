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

logger = logging.getLogger(__name__)

def call_subparser(parser):
    """
    Create a subparser for calling.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("xg_path", type=str,
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
    parser.add_argument("--gams", nargs='+', required=True,
                        help="GAMs to call.  One per chromosome. Must be same length/order as --chroms")

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
    parser.add_argument("--genotype_opts", type=str,
                        help="arguments to pass to vg genotype (wrapped in \"\")")
    parser.add_argument("--filter_opts", type=str,
                        help="argument to pass to vg filter (wrapped in \"\")")
    parser.add_argument("--calling_cores", type=int,
                        help="number of threads during the variant calling step")
    parser.add_argument("--vcf_offsets", nargs='+', default=[],
                         help="offset(s) to apply to output vcfs(s). (order of --chroms)")

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

def run_vg_call(job, options, vg_id, gam_id, pileup_id = None, xg_id = None,
                path_names = [], seq_names = [], seq_offsets = [], seq_lengths = [],
                filter_opts = [], pu_opts = [], call_opts = [],
                keep_pileup = False, keep_xg = False, keep_gam = False,
                keep_augmented = False, name = 'call'):
    """ Run vg call on a single graph.

    Returns (vcf_id, pileup_id, xg_id, gam_id, augmented_graph_id).  pileup_id and xg_id
    can be same as input if they are not computed.  If pileup/xg/augmented are 
    computed, the returned ids will be None unless appropriate keep_flag set
    (to prevent sending them to the file store if they aren't wanted)

    User is responsible to make sure that options passed in options.*_opts don't conflict
    with seq_names, seq_offsets, seq_lengths etc. If not provided, the pileup is computed.

    gam filtering is only done if filter_opts are passed in. 

    name option is only for working filenames (to make more readable)

    """

    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    vg_path = os.path.join(work_dir, '{}.vg'.format(name))
    read_from_store(job, options, vg_id, vg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(name))
    read_from_store(job, options, gam_id, gam_path)
    xg_path = os.path.join(work_dir, '{}.xg'.format(name))
    defray = filter_opts and ('-D' in filter_opts or '--defray-ends' in filter_opts)
    if xg_id and defray:
        read_from_store(job, options, xg_id, xg_path)
    pu_path = os.path.join(work_dir, '{}.pu'.format(name))
    if pileup_id:
        read_from_store(job, options, pileup_id, pu_path)
        
    # we only need an xg if using vg filter -D
    if not xg_id and defray:
        options.drunner.call(job, ['vg', 'index', os.path.basename(vg_path), '-x',
                                   os.path.basename(xg_path), '-t', str(options.calling_cores)],
                             work_dir = work_dir)
        if keep_xg:
            xg_id = write_to_store(job, options, xg_path)

    # pileup and filter together
    if filter_opts or not pileup_id:
        gam_filter_path = gam_path + '.filter'            
        if filter_opts:
            command = [['vg', 'filter', os.path.basename(gam_path), '-t', '1'] + filter_opts]
            if defray:
                command[0] += ['-x', os.path.basename(xg_path)]
            if keep_gam:
                command.append(['tee', os.path.basename(gam_filter_path)])
        else:
            command = [['cat', os.path.basename(gam_path)]]
        command.append(['vg', 'pileup', os.path.basename(vg_path), '-',
                        '-t', str(options.calling_cores)] + pu_opts)
        with open(pu_path, 'w') as pu_stream:
            options.drunner.call(job, command, work_dir=work_dir, outfile=pu_stream)
        if keep_gam and filter_opts:
            gam_id = write_to_store(job, options, gam_filter_path)
        if keep_pileup:
            pileup_id = write_to_store(job, options, pu_path)
        
    # call
    try:
        vcf_path = os.path.join(work_dir, '{}_call.vcf'.format(name))
        vcf_log_path = os.path.join(work_dir, '{}_call_log.txt'.format(name))
        aug_graph_id = None
        
        with open(vcf_path, 'w') as vgcall_stdout, open(vcf_log_path, 'w') as vgcall_stderr:
            command = ['vg', 'call', os.path.basename(vg_path), os.path.basename(pu_path), '-t',
                     str(options.calling_cores)] + ['-S', options.sample_name]
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
            aug_path = os.path.join(work_dir, '{}_aug.vg'.format(name))
            if keep_augmented:
                command.append(['-A', os.path.basename(aug_path)])
            options.drunner.call(job, command, work_dir=work_dir,
                                 outfile=vgcall_stdout, errfile=vgcall_stderr)
            if keep_augmented:
                aug_graph_id = write_to_store(job, options, aug_path)

        vcf_id = write_to_store(job, options, vcf_path)

    except Exception as e:
        logging.error("Failed. Dumping files.")
        write_to_store(job, options, vg_path, True)
        write_to_store(job, options, pu_path, True)
        write_to_store(job, options, vcf_path + ".call_log.txt", True)
        raise e
        
    return vcf_id, pileup_id, xg_id, gam_id, aug_graph_id


def run_vg_genotype(job, options, vg_id, gam_id, xg_id = None,
                    path_names = [], seq_names = [], seq_offsets = [],
                    seq_lengths = [], genotype_opts = [], filter_opts = [],
                    keep_xg = False, keep_gam = False, keep_augmented = False,
                    name = 'genotype'):
    """ Run vg genotype on a single graph.

    Returns (vcf_id, xg_id, gam_id, augmented_graph_id).  xg_id
    can be same as input if not computed.  If xg/augmented are 
    computed, the returned ids will be None unless appropriate keep_flag set
    (to prevent sending them to the file store if they aren't wanted)

    User is responsible to make sure that options passed in options.*_opts don't conflict
    with seq_names, seq_offsets, seq_lengths etc. If not provided, the pileup is computed.

    gam filtering is only done if filter_opts are passed in. 

    name option is only for working filenames (to make more readable)
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    vg_path = os.path.join(work_dir, '{}.vg'.format(name))
    read_from_store(job, options, vg_id, vg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(name))
    read_from_store(job, options, gam_id, gam_path)
    xg_path = os.path.join(work_dir, '{}.xg'.format(name))
    defray = filter_opts and ('-D' in filter_opts or '--defray-ends' in filter_opts)
    if xg_id and defray:
        read_from_store(job, options, xg_id, xg_path)
        
    # we only need an xg if using vg filter -D
    if not xg_id and defray:
        options.drunner.call(job, ['vg', 'index', os.path.basename(vg_path), '-x',
                                   os.path.basename(xg_path), '-t', str(options.calling_cores)],
                             work_dir = work_dir)
        if keep_xg:
            xg_id = write_to_store(job, options, xg_path)

    # filter
    if filter_opts:
        gam_filter_path = gam_path + '.filter'
        with open(gam_filter_path, 'w') as gam_filter_stream:
            command = [['vg', 'filter', os.path.basename(gam_path), '-t', '1'] + filter_opts]
            if defray:
                command[0] += ['-x', os.path.basename(xg_path)]
            options.drunner.call(job, command, work_dir=work_dir, outfile=gam_filter_stream)
        if keep_gam:
            gam_id = write_to_store(job, options, gam_filter_path)
    else:
        gam_filter_path = gam_path    

    # index the gam
    gam_index_path = os.path.join(work_dir, gam_filter_path + '.index')
    command = ['vg', 'index', '-N', os.path.basename(gam_filter_path),
               '-d', os.path.basename(gam_index_path)]
    options.drunner.call(job, command, work_dir=work_dir)

    # genotype
    try:
        vcf_path = os.path.join(work_dir, '{}_genotype.vcf'.format(name))
        vcf_log_path = os.path.join(work_dir, '{}_genotype_log.txt'.format(name))
        aug_graph_id = None
        
        with open(vcf_path, 'w') as vgcall_stdout, open(vcf_log_path, 'w') as vgcall_stderr:
            command = ['vg', 'genotype', os.path.basename(vg_path),
                       os.path.basename(gam_index_path), '-t',
                       str(options.calling_cores)] + genotype_opts + ['-s', options.sample_name, '-v']
            for path_name in path_names:
                command += ['-r', path_name]
            for seq_name in seq_names:
                command += ['-c', seq_name]
            for seq_length in seq_lengths:
                command += ['-l', seq_length]
            for seq_offset in seq_offsets:
                command += ['-o', seq_offset]
            aug_path = os.path.join(work_dir, '{}_aug.vg'.format(name))
            if keep_augmented:
                command.append(['-a', os.path.basename(aug_path)])
            options.drunner.call(job, command, work_dir=work_dir,
                                 outfile=vgcall_stdout, errfile=vgcall_stderr)
            if keep_augmented:
                aug_graph_id = write_to_store(job, options, aug_path)

        vcf_id = write_to_store(job, options, vcf_path)

    except Exception as e:
        logging.error("Failed. Dumping files.")
        write_to_store(job, options, vg_path, True)
        write_to_store(job, options, vcf_path + ".genotype_log.txt", True)
        raise e
        
    return vcf_id, xg_id, gam_id, aug_graph_id
        

def call_chunk(job, options, path_name, chunk_i, num_chunks, chunk_offset, clipped_chunk_offset,
               xg_file_id, vg_chunk_file_id, gam_chunk_file_id, path_size, vcf_offset):
    """ create VCF from a given chunk """
   
    RealtimeLogger.info("Running call_chunk on path {} and chunk {}".format(path_name, chunk_i))
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # output vcf path
    vcf_path = os.path.join(work_dir, 'chunk_{}_{}.vcf'.format(path_name, chunk_offset))

    # Run vg call
    if options.genotype:
        vcf_id, xg_id, gam_id, aug_graph_id = run_vg_genotype(
            job, options, vg_chunk_file_id, gam_chunk_file_id, xg_id=xg_file_id,
            path_names = [path_name],
            seq_names = [path_name],
            seq_offsets = [chunk_offset + vcf_offset],
            seq_lengths = [path_size],
            filter_opts = options.filter_opts,
            genotype_opts = options.genotype_opts,
            name = 'chunk_{}_{}'.format(path_name, chunk_offset))
    else:
        vcf_id, pu_id, xg_id, gam_id, aug_graph_id = run_vg_call(
            job, options, vg_chunk_file_id, gam_chunk_file_id, pileup_id = None, xg_id = xg_file_id,
            path_names = [path_name], 
            seq_names = [path_name],
            seq_offsets = [chunk_offset + vcf_offset],
            seq_lengths = [path_size],
            filter_opts = options.filter_opts, pu_opts = options.pileup_opts,
            call_opts = options.call_opts,
            name = 'chunk_{}_{}'.format(path_name, chunk_offset))
        
    # is this an unecessary copy or does toil cache handle it?
    read_from_store(job, options, vcf_id, vcf_path + '.us')

    # Sort the output
    sort_vcf(job, options.drunner, vcf_path + '.us', vcf_path)
    command=['bgzip', '{}'.format(os.path.basename(vcf_path))]
    options.drunner.call(job, command, work_dir=work_dir)
    command=['tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
    options.drunner.call(job, command, work_dir=work_dir)
    
    # do the vcf clip
    left_clip = 0 if chunk_i == 0 else options.overlap / 2
    right_clip = 0 if chunk_i == num_chunks - 1 else options.overlap / 2
    clip_path = os.path.join(work_dir, 'chunk_{}_{}_clip.vcf'.format(path_name, chunk_offset))
    with open(clip_path, "w") as clip_path_stream:
        offset = vcf_offset + 1
        command=['bcftools', 'view', '-r', '{}:{}-{}'.format(
            path_name, offset + clipped_chunk_offset + left_clip,
            offset + clipped_chunk_offset + options.call_chunk_size - right_clip - 1),
                 os.path.basename(vcf_path) + ".gz"]
        options.drunner.call(job, command, work_dir=work_dir, outfile=clip_path_stream)

    # save clip.vcf files to job store
    clip_file_id = write_to_store(job, options, clip_path)
    
    return clip_file_id

def run_all_calling(job, options, xg_file_id, chr_gam_ids, chroms):
    """
    Call all the chromosomes and return a merged up vcf/tbi pair
    """
    vcf_tbi_file_id_pair_list = []
    assert len(chr_gam_ids) > 0
    for i in range(len(chr_gam_ids)):
        alignment_file_id = chr_gam_ids[i]
        if len(chr_gam_ids) > 1:
            # 1 gam per chromosome
            chr_label = [chroms[i]]
            chr_offset = [options.vcf_offsets[i]] if options.vcf_offsets else [0]
        else:
            # single gam with one or more chromosomes
            chr_label = chroms
            chr_offset = options.vcf_offsets if options.vcf_offsets else [0] * len(chroms)
        vcf_tbi_file_id_pair = job.addChildJobFn(run_calling, options, xg_file_id,
                                                 alignment_file_id, chr_label, chr_offset,
                                                 cores=options.call_chunk_cores,
                                                 memory=options.call_chunk_mem,
                                                 disk=options.call_chunk_disk).rv()
        vcf_tbi_file_id_pair_list.append(vcf_tbi_file_id_pair)

    return job.addFollowOnJobFn(run_merge_vcf, options, vcf_tbi_file_id_pair_list,
                                cores=options.call_chunk_cores,
                                memory=options.call_chunk_mem,
                                disk=options.call_chunk_disk).rv()

def run_merge_vcf(job, options, vcf_tbi_file_id_pair_list):
    """ Merge up a bunch of chromosome VCFs """

    RealtimeLogger.info("Completed gam merging and gam path variant calling.")
    RealtimeLogger.info("Starting vcf merging vcf files.")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    out_store = IOStore.get(options.out_store)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vcf_merging_file_key_list = [] 
    for i, vcf_tbi_file_id_pair in enumerate(vcf_tbi_file_id_pair_list):
        vcf_file = os.path.join(work_dir, 'vcf_chunk_{}.vcf.gz'.format(i))
        vcf_file_idx = '{}.tbi'.format(vcf_file)
        read_from_store(job, options, vcf_tbi_file_id_pair[0], vcf_file)
        read_from_store(job, options, vcf_tbi_file_id_pair[1], vcf_file_idx)
        vcf_merging_file_key_list.append(os.path.basename(vcf_file))

    vcf_merged_file_key = "" 
    if len(vcf_merging_file_key_list) > 1:
        # merge vcf files
        vcf_merged_file_key = "{}.vcf.gz".format(options.sample_name)
        command = ['bcftools', 'concat', '-O', 'z', '-o', os.path.basename(vcf_merged_file_key)]
        command +=  vcf_merging_file_key_list
        options.drunner.call(job, command, work_dir=work_dir)
        command=['bcftools', 'tabix', '-f', '-p', 'vcf', os.path.basename(vcf_merged_file_key)]
        options.drunner.call(job, command, work_dir=work_dir)
    else:
        vcf_merged_file_key = vcf_merging_file_key_list[0]

    # save variant calling results to the output store
    out_store_key = "{}.vcf.gz".format(options.sample_name)
    vcf_file = os.path.join(work_dir, vcf_merged_file_key)
    vcf_file_idx = vcf_file + ".tbi"
    
    vcf_file_id = write_to_store(job, options, vcf_file)
    vcf_idx_file_id = write_to_store(job, options, vcf_file_idx)
    
    # checkpoint to output store if not done above
    if not options.force_outstore:
        write_to_store(job, options, vcf_file, use_out_store = True, out_store_key = out_store_key)
        write_to_store(job, options, vcf_file_idx, use_out_store = True, out_store_key = out_store_key + ".tbi")

    return vcf_file_id, vcf_idx_file_id


def run_calling(job, options, xg_file_id, alignment_file_id, path_names, vcf_offsets):
    
    RealtimeLogger.info("Running variant calling on path(s) {} from alignment file {}".format(','.join(path_names), str(alignment_file_id)))
        
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Tame for work files
    tag = path_names[0] if len(path_names) == 1 else 'chroms'

    # Download the input from the store
    xg_path = os.path.join(work_dir, 'graph.vg.xg')
    read_from_store(job, options, xg_file_id, xg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(options.sample_name))
    read_from_store(job, options, alignment_file_id, gam_path)

    # Index the gam
    gam_index_path = gam_path + '.index'
    index_cmd = ['vg', 'index', '-N', os.path.basename(gam_path), '-d',
                 os.path.basename(gam_index_path), '-t', str(options.gam_index_cores)]
    options.drunner.call(job, index_cmd, work_dir=work_dir)

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
                 '-a', os.path.basename(gam_index_path), '-c', str(options.chunk_context),
                 '-P', os.path.basename(path_list),
                 '-g',
                 '-s', str(options.call_chunk_size),
                 '-o', str(options.overlap),
                 '-b', 'call_chunk_{}'.format(tag),
                 '-t', str(options.call_chunk_cores),
                 '-R', os.path.basename(output_bed_chunks_path)]
    options.drunner.call(job, chunk_cmd, work_dir=work_dir)

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
        
    # Go through the BED output of vg chunk, adding a child calling job for
    # each chunk                
    clip_file_ids = []
    for toks in bed_lines:
        chunk_bed_chrom = toks[0]
        chunk_bed_start = int(toks[1])
        chunk_bed_end = int(toks[2])
        chunk_bed_size = chunk_bed_end - chunk_bed_start
        gam_chunk_path = os.path.join(work_dir, os.path.basename(toks[3].strip()))
        assert os.path.splitext(gam_chunk_path)[1] == '.gam'
        vg_chunk_path = os.path.splitext(gam_chunk_path)[0] + '.vg'
        gam_chunk_file_id = write_to_store(job, options, gam_chunk_path)
        vg_chunk_file_id = write_to_store(job, options, vg_chunk_path)
        chunk_i = cur_path_offset[chunk_bed_chrom]        
        clipped_chunk_offset = chunk_i * options.call_chunk_size - chunk_i * options.overlap
        cur_path_offset[chunk_bed_chrom] += 1
        
        clip_file_id = job.addChildJobFn(call_chunk, options, chunk_bed_chrom, chunk_i,
                                         len(bed_lines),
                                         chunk_bed_start, clipped_chunk_offset,
                                         xg_file_id, vg_chunk_file_id, gam_chunk_file_id,
                                         path_size[chunk_bed_chrom], offset_map[chunk_bed_chrom],
                                         cores=options.calling_cores,
                                         memory=options.calling_mem, disk=options.calling_disk).rv()
        clip_file_ids.append(clip_file_id)


    vcf_gz_tbi_file_id_pair = job.addFollowOnJobFn(merge_vcf_chunks, options, tag,
                                                   clip_file_ids,
                                                   cores=options.call_chunk_cores,
                                                   memory=options.call_chunk_mem,
                                                   disk=options.call_chunk_disk).rv()
 
    return vcf_gz_tbi_file_id_pair


def merge_vcf_chunks(job, options, path_name, clip_file_ids):
    """ merge a bunch of clipped vcfs created above, taking care to 
    fix up the headers.  everything expected to be sorted already """
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vcf_path = os.path.join(work_dir, path_name + ".vcf")
    
    for chunk_i, clip_file_id in enumerate(clip_file_ids):
        
        # Download clip.vcf file from the store
        clip_path = os.path.join(work_dir, 'clip_{}.vcf'.format(chunk_i))
        read_from_store(job, options, clip_file_id, clip_path)

        if chunk_i == 0:
            # copy everything including the header
            with open(vcf_path, "w") as outfile:
                options.drunner.call(job, ['cat', os.path.basename(clip_path)], outfile=outfile,
                                     work_dir=work_dir)
        else:
            # add on everythin but header
            with open(vcf_path, "a") as outfile:
                options.drunner.call(job, ['bcftools', 'view', '-H', os.path.basename(clip_path)],
                                     outfile=outfile, work_dir=work_dir)

    # add a compressed indexed version
    vcf_gz_file = vcf_path + ".gz"
    with open(vcf_gz_file, "w") as vcf_gz_file_stream:
        command=['bgzip', '-c', '{}'.format(os.path.basename(vcf_path))]
        options.drunner.call(job, command, work_dir=work_dir, outfile=vcf_gz_file_stream)
    command=['bcftools', 'tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
    options.drunner.call(job, command, work_dir=work_dir)

    # Save merged vcf files to the job store
    vcf_gz_file_id = write_to_store(job, options, vcf_path+".gz")
    vcf_tbi_file_id = write_to_store(job, options, vcf_path+".gz.tbi")

    RealtimeLogger.info("Completed variant calling on path {}".format(path_name))

    return vcf_gz_file_id, vcf_tbi_file_id
    
def call_main(options):
    """ no harm in preserving command line access to chunked_call for debugging """
    
    # make the docker runner
    options.drunner = ContainerRunner(
        container_tool_map = get_container_tool_map(options))

    require(len(options.chroms) == len(options.gams) or len(options.gams) == 1,
            'Number of --chroms must be 1 or same as number of --gams')
    require(not options.vcf_offsets or len(options.vcf_offsets) == len(options.chroms),
            'Number of --vcf_offsets if specified must be same as number of --chroms')
    
    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'call'

    # Throw error if something wrong with IOStore string
    IOStore.get(options.out_store)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()

            # Upload local files to the remote IO Store
            inputXGFileID = import_to_store(toil, options, options.xg_path)
            inputGamFileIDs = []
            for inputGamFileID in options.gams:
                inputGamFileIDs.append(import_to_store(toil, options, inputGamFileID))

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_all_calling, options, inputXGFileID, inputGamFileIDs,
                                     options.chroms,
                                     cores=options.misc_cores,
                                     memory=options.misc_mem,
                                     disk=options.misc_disk)

            
            # Run the job and store the returned list of output files to download
            toil.start(root_job)
        else:
            toil.restart()
                
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
    
