#!/usr/bin/env python2.7
"""
Generate a VCF from a GAM and XG by splitting into GAM/VG chunks.
"""

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
from toil_vg.vg_chunk import *
from toil_vg.vg_augment import *

logger = logging.getLogger(__name__)

def call_subparser(parser):
    """
    Create a subparser for calling.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified "
                        "using same syntax as toil jobStore")
    parser.add_argument("--gam", type=make_url, required=True,
                        help="GAM to call")
    parser.add_argument("--graph", type=make_url, required=True,
                        help="graph to augment")
    parser.add_argument("--snarls", type=make_url,
                        help="Path to snarls file")
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common calling options shared with toil_vg pipeline
    call_parse_args(parser)

    # Add common chunking options
    chunk_parse_args(parser, path_components=False)

    # Add common augmenting options
    augment_parse_args(parser)

    # Add common docker options shared with toil_vg pipeline
    add_container_tool_parse_args(parser)
                        

def call_parse_args(parser):
    """ centralize calling parameters here """
    parser.add_argument("--sample", type=str,
                        help="sample name for output VCF")
    parser.add_argument("--genotype_vcf", type=make_url,
                        help="genotype the given VCF.  Input graph must contain alt paths from VCF")
    parser.add_argument("--recall", action="store_true",
                        help="do not augment. only recall variants in graph. (on by default with --genotype_vcf)")
    parser.add_argument("--filter_opts", type=str,
                        help="argument to pass to vg filter (wrapped in \"\")")
    parser.add_argument("--vcf_offsets", nargs='+', default=[],
                        help="offset(s) to apply to output vcfs(s). (order of --chroms)")
    parser.add_argument("--min_call_support", type=float,
                        help="minimum support to make a call")
    parser.add_argument("--calling_cores", type=int,
                        help="number of threads during the variant calling step")
    parser.add_argument("--calling_mem", type=str,
                        help="memory alotment during the variant calling step")
    parser.add_argument("--gam_chunking", action="store_true",
                        help="split the GAM into chromosome components (in addition to the graph)")

def validate_call_options(options):
    require(not options.vcf_offsets or len(options.vcf_offsets) == len(options.ref_paths),
            'Number of --vcf_offsets if specified must be same as number of --ref_paths')
    require(not options.snarls or options.recall or options.genotype_vcf,
            '--snarls can only be used with --recall or --genotype_vcf')
    require(not options.connected_component_chunking or not options.genotype_vcf,
            '--only path chunking (not --connected_component_chunking) supported with --genotype_vcf')
    require(not options.connected_component_chunking or not options.ref_path_chunking,
            '--connected_components cannot be used with --ref_path_chunking')

def run_chunked_calling(job, context,
                        graph_id,
                        graph_basename,
                        gam_id,
                        gam_basename,
                        batch_input=None,
                        snarls_id=None,
                        genotype_vcf_id=None,
                        genotype_tbi_id=None,
                        sample=None,
                        augment=False,
                        connected_component_chunking=False,
                        output_format=None,
                        min_augment_coverage=None,
                        expected_coverage=None,
                        min_mapq=None,
                        min_baseq=None,
                        ref_paths=[],
                        ref_path_chunking=True,
                        min_call_support=None,
                        vcf_offsets={},
                        gam_chunking=False):

    # simple way to keep follow-ons down the tree
    child_job = Job()
    job.addChild(child_job)

    out_vcf_name = remove_ext(graph_basename)
    if sample:
        out_vcf_name += '_' + sample

    # base case: only one input
    if batch_input is None:
        # chunk if necessary
        if connected_component_chunking or ref_path_chunking:

            chunk_job = child_job.addChildJobFn(run_chunking, context,
                                                graph_id=graph_id,
                                                graph_basename=graph_basename,
                                                chunk_paths=ref_paths,
                                                connected_component_chunking=connected_component_chunking,
                                                output_format=output_format,
                                                gam_id=gam_id if gam_chunking else None,
                                                to_outstore=False,
                                                cores=context.config.chunk_cores,
                                                memory=context.config.chunk_mem,
                                                disk=context.config.chunk_disk)

            batch_input = chunk_job.rv()

            # recurse on chunks
            recurse_job = child_job.addFollowOnJobFn(run_chunked_calling, context,
                                                     graph_id=None,
                                                     graph_basename=graph_basename,
                                                     gam_id=gam_id,
                                                     gam_basename=gam_basename,
                                                     batch_input=batch_input,
                                                     snarls_id=snarls_id,
                                                     genotype_vcf_id=genotype_vcf_id,
                                                     genotype_tbi_id=genotype_tbi_id,
                                                     sample=sample,
                                                     augment=augment,
                                                     connected_component_chunking=connected_component_chunking,
                                                     output_format=output_format,
                                                     min_augment_coverage=min_augment_coverage,
                                                     expected_coverage=expected_coverage,
                                                     min_mapq=min_mapq,
                                                     min_baseq=min_baseq,
                                                     ref_paths=ref_paths,
                                                     ref_path_chunking=ref_path_chunking,
                                                     min_call_support=min_call_support,
                                                     vcf_offsets=vcf_offsets,
                                                     gam_chunking=gam_chunking)
            return recurse_job.rv()
        else:
            # convert if we're augmenting and not chunking
            if augment and os.path.splitext(graph_basename)[1] != '.' + output_format:
                convert_job = child_job.addChildJobFn(run_convert, context,
                                                      graph_id=graph_id,
                                                      graph_basename=graph_basename,
                                                      output_format=output_format,
                                                      disk=context.config.calling_disk)
                graph_id = convert_job.rv()
                graph_basename = os.path.splitext(graph_basename)[0] + '.' + output_format
                # todo: clean up
                next_job = Job()
                child_job.addFollowOn(next_job)
                child_job = next_job
                
            #phony up chunk output for single input
            batch_input = { 'all' : [graph_id, graph_basename] }
            if gam_id:
                batch_input['all'] += [gam_id, gam_basename]

    # run the calling on each chunk
    assert batch_input

    call_results = []
    in_gam_id = gam_id
    in_gam_basename = gam_basename
    for chunk_name, chunk_results in list(batch_input.items()):
        calling_root_job = Job()
        child_job.addChild(calling_root_job)

        graph_id = chunk_results[0]
        graph_basename = chunk_results[1]
        if gam_chunking:
            gam_id = chunk_results[2]
            gam_basename = chunk_results[3]
        else:
            gam_id = in_gam_id
            gam_basename = in_gam_basename
            
        if augment:
            augment_job = calling_root_job.addChildJobFn(run_augmenting, context,
                                                         graph_id=graph_id,
                                                         graph_basename=graph_basename,
                                                         gam_id=gam_id,
                                                         gam_basename=gam_basename,
                                                         augment_gam=True,
                                                         min_augment_coverage=min_augment_coverage,
                                                         expected_coverage=expected_coverage,
                                                         min_mapq=min_mapq,
                                                         min_baseq=min_baseq,
                                                         to_outstore=True,
                                                         cores=context.config.augment_cores,
                                                         memory=context.config.augment_mem,
                                                         disk=context.config.augment_disk)
            graph_id = augment_job.rv(0)
            graph_basename = os.path.splitext(graph_basename)[0] + '-aug' + os.path.splitext(graph_basename)[1]
            gam_id = augment_job.rv(1)
            gam_basename = os.path.splitext(gam_basename)[0] + '-aug' + os.path.splitext(gam_basename)[1]

        # When path chunking, we subset our reference paths down to the current path
        if ref_path_chunking:
            ref_path = [chunk_name]
        else:
            ref_path = ref_paths

        calling_job = calling_root_job.addFollowOnJobFn(run_calling, context,
                                                        graph_id=graph_id,
                                                        graph_basename=graph_basename,
                                                        gam_id=gam_id,
                                                        gam_basename=gam_basename,
                                                        snarls_id=snarls_id,
                                                        genotype_vcf_id=genotype_vcf_id,
                                                        genotype_tbi_id=genotype_tbi_id,
                                                        sample=sample,
                                                        expected_coverage=expected_coverage,
                                                        min_mapq=min_mapq,
                                                        ref_paths=ref_path,
                                                        min_call_support=min_call_support,
                                                        vcf_offsets=vcf_offsets,
                                                        to_outstore=False,
                                                        cores=context.config.calling_cores,
                                                        memory=context.config.calling_mem,
                                                        disk=context.config.calling_disk)

        call_results.append((chunk_name, calling_job.rv()))

    concat_job = child_job.addFollowOnJobFn(run_concat_vcfs, context,
                                            out_name = out_vcf_name,
                                            vcf_ids = None,
                                            tbi_ids = None,
                                            write_to_outstore = True,
                                            call_timers_lists = [],
                                            batch_data = call_results)

    return concat_job.rv()

def run_calling(job, context,
                graph_id,
                graph_basename,
                gam_id,
                gam_basename,
                snarls_id=None,
                genotype_vcf_id=None,
                genotype_tbi_id=None,
                sample=None,
                expected_coverage=None,
                min_mapq=None,
                ref_paths=None,
                min_call_support=None,
                vcf_offsets=None,
                to_outstore=True):
    """
    run vg pack and vg call to make a vcf
    """
    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    graph_path = os.path.join(work_dir, graph_basename)
    job.fileStore.readGlobalFile(graph_id, graph_path)    
    gam_path = os.path.join(work_dir, gam_basename)
    job.fileStore.readGlobalFile(gam_id, gam_path)
    snarls_path = remove_ext(graph_path) + '.snarls'
    if snarls_id:
        job.fileStore.readGlobalFile(snarls_id, snarls_path)
    genotype_vcf_path = remove_ext(graph_path) + '.vcf.gz'
    genotype_tbi_path = genotype_vcf_path + '.tbi'
    if genotype_vcf_id:
        job.fileStore.readGlobalFile(genotype_vcf_id, genotype_vcf_path)
        job.fileStore.readGlobalFile(genotype_tbi_id, genotype_tbi_path)
        # (I think this is the only case where we want to do this)
        if ref_paths and len(ref_paths) == 1:
            # cut the VCF to match our chunk
            sub_vcf_path = remove_ext(genotype_vcf_path, '.vcf.gz') + '_' + ref_paths[0] + '.vcf.gz'
            sub_tbi_path = sub_vcf_path + '.tbi'
            with open(sub_vcf_path, 'wb') as sub_vcf_file:
                context.runner.call(job, ['bcftools', 'view', os.path.basename(genotype_vcf_path), '-r',
                                          ref_paths[0], '-O', 'z'], work_dir = work_dir, outfile = sub_vcf_file)
            context.runner.call(job, ['tabix', '-f', '-p', 'vcf', os.path.basename(sub_vcf_path)], work_dir = work_dir)
            genotype_vcf_path, genotype_tbi_path = sub_vcf_path, sub_tbi_path
            
    # Run vg pack
    pack_path = remove_ext(gam_path) + '.pack'

    pack_cmd = ['vg', 'pack', '-x', os.path.basename(graph_path), '-g', os.path.basename(gam_path),
                '-o', os.path.basename(pack_path), '-t', str(job.cores)]

    if min_mapq is not None:
        pack_cmd += ['-Q', str(min_mapq)]
    if expected_coverage is not None:
        pack_cmd += ['-c', str(expected_coverage)]

    try:
        context.runner.call(job, pack_cmd, work_dir = work_dir)
    except Exception as e:
        logging.error("Pack failed. Dumping input files to outstore.")
        for dump_path in [graph_path, gam_path]:
            if dump_path and os.path.isfile(dump_path):
                context.write_output_file(job, dump_path)
        raise
    
    # Run vg call
    out_vcf_path = remove_ext(graph_path)
    if sample:
        out_vcf_path += '_' + sample
    out_vcf_path += '.vcf.gz'

    call_cmd = ['vg', 'call', os.path.basename(graph_path), '-k', os.path.basename(pack_path), '-t', str(job.cores)]
    if genotype_vcf_id:
        call_cmd += ['-v', os.path.basename(genotype_vcf_path)]
    if snarls_id:
        call_cmd += ['-r', os.path.basename(snarls_path)]
    if ref_paths:
        for ref_path in ref_paths:
            call_cmd += ['-p', ref_path]
            if vcf_offsets and ref_path in vcf_offsets:
                call_cmd += ['-o', str(vcf_offsets[ref_path])]
                
    if sample:
        call_cmd += ['-s', sample]
    if context.config.call_opts:
        call_cmd += context.config.call_opts

    try:
        with open(out_vcf_path, 'wb') as out_vcf_file:
            context.runner.call(job, [call_cmd, ['bgzip']], work_dir = work_dir, outfile = out_vcf_file)
    except Exception as e:
        logging.error("Pack failed. Dumping input files to outstore.")
        for dump_path in [graph_path, pack_path, snarls_path, genotype_vcf_path, genotype_tbi_path]:
            if dump_path and os.path.isfile(dump_path):
                context.write_output_file(job, dump_path)
        raise
        
    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', os.path.basename(out_vcf_path)], work_dir=work_dir)

    write_fn = context.write_output_file if to_outstore else context.write_intermediate_file
    return write_fn(job, out_vcf_path), write_fn(job, out_vcf_path + '.tbi')

def run_convert(job, context,
                graph_id,
                graph_basename,
                output_format):
    """ convert a graph """
    
    work_dir = job.fileStore.getLocalTempDir()

    # Read our input files from the store
    graph_path = os.path.join(work_dir, graph_basename)
    job.fileStore.readGlobalFile(graph_id, graph_path)

    out_graph_path = os.path.splitext(graph_path)[0] + '.' + output_format

    with open(out_graph_path, 'wb') as out_graph_file:
        convert_cmd = ['vg', 'convert', os.path.basename(graph_path)]
        flag_map = { 'vg' : '-v', 'hg' : '-a', 'pg' : '-p' }
        assert output_format in flag_map
        convert_cmd += [flag_map[output_format]]

        context.runner.call(job, convert_cmd, work_dir=work_dir, outfile = out_graph_file)

    return context.write_intermediate_file(job, out_graph_path)

def run_concat_vcfs(job, context, out_name, vcf_ids, tbi_ids = [], write_to_outstore = False,
                    call_timers_lists = [], batch_data=None):
    """ Concat up a bunch of VCFs. Input assumed to be bgzipped iff tbi_ids specified """

    # hack to re-use this function with new interface
    if batch_data:
        vcf_ids = [result[1][0] for result in batch_data]
        tbi_ids = [result[1][1] for result in batch_data]
        
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

def run_filtering(job, context,
                  graph_id,
                  graph_basename,
                  gam_id, 
                  gam_basename,
                  filter_opts):

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Don't do anything if there's nothing to pass to filter
    if not filter_opts:
        return gam_id

    if not graph_basename:
        graph_basename = 'graph'
    if not gam_basename:
        gam_basename = 'aln.gam'
    
    # Download the input
    gam_path = os.path.join(work_dir, gam_basename)
    job.fileStore.readGlobalFile(gam_id, gam_path)

    if '-D' in filter_opts or '--defray-ends' in filter_opts:
        graph_path = os.path.join(work_dir, 'graph')
        job.fileStore.readGlobalFile(graph_id, graph_path)
        filter_opts += ['-x', os.path.basename(graph_path)]

    filter_path = os.path.splitext(gam_path)[0] + '-filter.gam'
    with open(filter_path, 'wb') as filter_file:
        context.runner.call(job, ['vg', 'filter', os.path.basename(gam_path)] + filter_opts,
                            work_dir=work_dir, outfile = filter_file)

    return context.write_intermediate_file(job, filter_path)
    
def call_main(context, options):
    """
    Wrapper for vg filter / pack / call
    """

    validate_call_options(options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)

            # Upload local files to the job store
            inputGraphFileID = importer.load(options.graph)
            inputGamFileID = importer.load(options.gam)
            snarlsFileID = None
            if options.snarls:
                snarlsFileID = importer.load(options.snarls)
            inputVcfID = None
            inputTbiID = None
            if options.genotype_vcf:
                inputVcfID = importer.load(options.genotype_vcf)
                inputTbiID = importer.load(options.genotype_vcf + '.tbi')

            importer.wait()

            # transform our vcf_offsets into a dict
            if options.vcf_offsets:
                vcf_offset_dict = {}
                for (ref_path_name, vcf_offset) in zip(options.ref_paths, options.vcf_offsets):
                    vcf_offset_dict[ref_path_name] = vcf_offset
                options.vcf_offsets = vcf_offset_dict

            # Filter the gam as a preprocessing step
            # It's arguable this would be more efficient, especially in terms of disk, if
            # it was moved after chunking, but it makes the logic more complicated and I want
            # to eventually get away from running it at all.
            root_job = None
            if context.config.filter_opts:
                root_job = Job.wrapJobFn(run_filtering, context,
                                         graph_id = importer.resolve(inputGraphFileID),
                                         graph_basename = os.path.basename(options.graph),
                                         gam_id = importer.resolve(inputGamFileID),
                                         gam_basename = os.path.basename(options.gam),
                                         filter_opts = context.config.filter_opts,
                                         cores=context.config.calling_cores,
                                         memory=context.config.calling_mem,
                                         disk=context.config.calling_disk)
                filtered_gam_id = root_job.rv()
                filtered_gam_basename = os.path.splitext(os.path.basename(options.gam))[0] + '.filter.gam'
            else:
                filtered_gam_id = importer.resolve(inputGamFileID)
                filtered_gam_basename = os.path.basename(options.gam)

            calling_job = Job.wrapJobFn(run_chunked_calling, context,
                                        graph_id = importer.resolve(inputGraphFileID),
                                        graph_basename = os.path.basename(options.graph),
                                        gam_id = filtered_gam_id,
                                        gam_basename = filtered_gam_basename,
                                        batch_input=None,
                                        snarls_id = importer.resolve(snarlsFileID),
                                        genotype_vcf_id = importer.resolve(inputVcfID),
                                        genotype_tbi_id = importer.resolve(inputTbiID),
                                        sample=options.sample,
                                        connected_component_chunking=options.connected_component_chunking,
                                        augment=not options.recall and options.genotype_vcf is None,
                                        output_format=options.output_format,
                                        min_augment_coverage=options.min_augment_coverage,
                                        expected_coverage=options.expected_coverage,
                                        min_mapq=options.min_mapq,
                                        min_baseq=options.min_baseq,
                                        ref_paths=options.ref_paths,
                                        ref_path_chunking=options.ref_path_chunking,
                                        min_call_support=options.min_call_support,
                                        vcf_offsets=options.vcf_offsets,
                                        gam_chunking=options.gam_chunking)

            if root_job:
                root_job.addFollowOn(calling_job)
            else:
                root_job = calling_job

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
    
    
