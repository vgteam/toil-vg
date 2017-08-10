#!/usr/bin/env python2.7
"""
vg_construct.py: construct a graph from a vcf and fasta

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import getpass
import pdb
import logging

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context

logger = logging.getLogger(__name__)

def construct_subparser(parser):
    """
    Create a subparser for indexing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")

    parser.add_argument("--fasta", required=True, type=make_url,
                        help="Reference sequence in fasta format")
    parser.add_argument("--vcf", default=None, type=make_url,
                        help="Variants to make graph from")
    parser.add_argument("--regions", nargs='+',
                        help="1-based inclusive VCF coordinates in the form SEQ:START-END")
    parser.add_argument("--max_node_size", type=int,
                        help="Maximum node length")
    parser.add_argument("--alt_paths", action="store_true",
                        help="Save paths for alts with variant ID (required for GPBWT)")
    parser.add_argument("--flat_alts", action="store_true",
                        help="flat alts")
    parser.add_argument("--control_sample", default=None,
                        help="Make a positive and negative control using this sample")
    parser.add_argument("--construct_cores", type=int,
                        help="Number of cores for vg construct")
    parser.add_argument("--graph_name",
                        help="Name of output graph.  If specified with multiple regions, they will be merged together")
    

    # Add common docker options
    add_container_tool_parse_args(parser)

    
def run_construct_with_controls(job, context, sample, fasta_id, fasta_name, vcf_id, vcf_name, tbi_id,
                                max_node_size, alt_paths, flat_alts, regions, sort_ids = True,
                                join_ids = True, merge_output_name=None):
    """ 
    construct a genome graph from a vcf.  also extract a negative and postive control using 
    a given sample and construct those too
    return ((vg ids), (positive control ids), (negative control ids))
    """
    if sample:
        make_controls = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id, sample)
        pos_control_vcf_id, pos_control_tbi_id = make_controls.rv(0), make_controls.rv(1)
        neg_control_vcf_id, neg_control_tbi_id = make_controls.rv(2), make_controls.rv(3)

        vcf_base = os.path.basename(vcf_name.rstrip('.gz').rstrip('.vcf'))
        pos_control_vcf_name = '{}_{}.vcf.gz'.format(vcf_base, sample)
        neg_control_vcf_name = '{}_minus_{}.vcf.gz'.format(vcf_base, sample)
        if regions:
            pos_region_names = [c.replace(':','-') + '_{}'.format(sample) for c in regions]
            neg_region_names = [c.replace(':','-') + '_minus_{}'.format(sample) for c in regions]
        else:
            pos_region_names = None
            neg_region_names = None
        if merge_output_name:
            pos_output_name = mergeoutput.name.rstrip('.vg') + '_{}.vg'.format(sample)
            neg_output_name = mergeoutput.name.rstrip('.vg') + '_minus_{}.vg'.format(sample)
        else:
            pos_output_name = None
            neg_output_name = None
        
        
        pos_control_vg_ids = job.addFollowOnJobFn(run_construct_genome_graph, context, fasta_id,
                                               fasta_name, pos_control_vcf_id, pos_control_vcf_name,
                                               pos_control_tbi_id,
                                               max_node_size, alt_paths, flat_alts, regions,
                                               pos_region_names,
                                               sort_ids=sort_ids, join_ids=join_ids,
                                               merge_output_name=pos_output_name).rv()

        neg_control_vg_ids = job.addFollowOnJobFn(run_construct_genome_graph, context, fasta_id,
                                               fasta_name, neg_control_vcf_id, neg_control_vcf_name,
                                               neg_control_tbi_id,
                                               max_node_size, alt_paths, flat_alts, regions,
                                               neg_region_names,
                                               sort_ids=sort_ids, join_ids=join_ids,
                                               merge_output_name=neg_output_name).rv()
    else:
        pos_control_vg_ids = None
        neg_control_vg_ids = None

    vg_ids = job.addChildJobFn(run_construct_genome_graph, context, fasta_id,
                               fasta_name, vcf_id, vcf_name, tbi_id,
                               max_node_size, alt_paths, flat_alts, regions,
                               [c.replace(':','-') for c in regions] if regions else None,
                               sort_ids=sort_ids, join_ids=join_ids,
                               merge_output_name=merge_output_name).rv()

    
    return (vg_ids, pos_control_vg_ids, neg_control_vg_ids)
                

def run_construct_genome_graph(job, context, fasta_id, fasta_name, vcf_id, vcf_name, tbi_id,
                              max_node_size, alt_paths, flat_alts, regions, region_names,
                              sort_ids = True, join_ids = True, merge_output_name=None):
    """ construct graph(s) from several regions in parallel.  we could eventually generalize this
    to accept multiple vcfs and/or fastas if needed, as well as to determine regions from file,
    but for now we only accept single files, and require region list.
    """

    work_dir = job.fileStore.getLocalTempDir()

    if not regions:
        regions, region_names = [None], ['genome']
        
    region_graph_ids = []    
    for region, region_name in zip(regions, region_names):
        region_graph_ids.append(job.addChildJobFn(run_construct_region_graph, context,
                                                  fasta_id, fasta_name,
                                                  vcf_id, vcf_name, tbi_id, region, region_name,
                                                  max_node_size, alt_paths, flat_alts,
                                                  sort_ids=sort_ids).rv())

    return job.addFollowOnJobFn(run_join_graphs, context, region_graph_ids, join_ids,
                                region_names, merge_output_name).rv()


def run_join_graphs(job, context, region_graph_ids, join_ids, region_names, merge_output_name = None):
    """
    join the ids of some graphs.  if a merge_output_name is given, cat them all together as well
    this function saves output to the outstore
    """
        
    work_dir = job.fileStore.getLocalTempDir()

    # download graph for each region
    region_files = []
    for region_graph_id, region_name in zip(region_graph_ids, region_names):
        region_file = '{}.vg'.format(region_name)
        job.fileStore.readGlobalFile(region_graph_id, os.path.join(work_dir, region_file), mutable=True)
        region_files.append(region_file)

    # if there's nothing to do, just write the files and return
    if len(region_graph_ids) == 1 or not join_ids:
        out_ids = []
        for region_file in region_files:
            out_ids.append(context.write_output_file(job, os.path.join(work_dir, region_file)))
        return out_ids

    if join_ids:
        # join the ids
        cmd = ['vg', 'ids', '-j'] + region_files
        context.runner.call(job, cmd, work_dir=work_dir)

    if merge_output_name is not None:
        if not merge_output_name.endswith('.vg'):
            merge_output_name += '.vg'
        assert merge_output_name[:-3] not in region_names
        with open(os.path.join(work_dir, merge_output_name), 'w') as merge_file:
            for region_file in region_files:
                with open(os.path.join(work_dir, region_file)) as cf:
                    shutil.copyfileobj(cf, merge_file)
        return [context.write_output_file(job, os.path.join(work_dir, merge_output_name))]
    else:
        return [context.write_output_file(job, os.path.join(work_dir, f)) for f in region_files]
        
    
def run_construct_region_graph(job, context, fasta_id, fasta_name, vcf_id, vcf_name, tbi_id,
                               region, region_name, max_node_size, alt_paths, flat_alts,
                               is_chrom = False, sort_ids = True):
    """ construct a graph from the vcf for a given region and return its id """

    work_dir = job.fileStore.getLocalTempDir()

    # Download input files
    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    job.fileStore.readGlobalFile(fasta_id, fasta_file)
    if vcf_id:
        vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
        job.fileStore.readGlobalFile(vcf_id, vcf_file)
        job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    RealtimeLogger.info("Fasta file {}".format(fasta_file))

    cmd = ['vg', 'construct', '-r', os.path.basename(fasta_file)]
    if vcf_id:
        cmd += ['-v', os.path.basename(vcf_file)]
    if region:
        cmd += ['-R', region]
    if is_chrom:
        cmd += ['-C']
    if max_node_size:
        cmd += ['-m', context.max_node_size]
    if alt_paths:
        cmd += '-a'
    if job.cores:
        cmd += ['-t', job.cores]

    if sort_ids:
        cmd = [cmd, ['vg', 'ids', '-s', '-']]

    vg_path = os.path.join(work_dir, region_name)
    with open(vg_path, 'w') as vg_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile = vg_file)

    return context.write_intermediate_file(job, vg_path)

def run_make_control_vcfs(job, context, vcf_id, vcf_name, tbi_id, sample):
    """ make a positive and negative control vcf 
    The positive control has only variants in the sample, the negative
    control has only variants not in the sample
    """

    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    # filter down to sample in question
    cmd = [['bcftools', 'view', os.path.basename(vcf_file), '-s', sample]]
    
    # remove anything that's not alt (probably cleaner way to do this0
    gfilter = 'GT="0" || GT="0|0" || GT="0/0"'
    gfilter += ' || GT="." || GT=".|." || GT="./."'
    gfilter += ' || GT=".|0" || GT="0/."'
    gfilter += ' || GT="0|." || GT="./0"'
    
    cmd.append(['bcftools', 'view', '-', '-O', 'z', '-e', gfilter])

    out_pos_name = os.path.basename(vcf_name)
    if out_pos_name.endswith('.gz'):
        out_pos_name = out_pos_name[:-3]
    if out_pos_name.endswith('.vcf'):
        out_pos_name = out_pos_name[:-4]
    out_neg_name = out_pos_name + '_minus_{}.vcf.gz'.format(sample)
    out_pos_name += '_{}.vcf.gz'.format(sample)

    with open(os.path.join(work_dir, out_pos_name), 'w') as out_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile = out_file)

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', out_pos_name], work_dir=work_dir)

    pos_control_vcf_id = context.write_output_file(job, os.path.join(work_dir, out_pos_name))
    pos_control_tbi_id = context.write_output_file(job, os.path.join(work_dir, out_pos_name + '.tbi'))

    # subtract the positive control to make the negative control
    cmd = ['bcftools', 'isec', os.path.basename(vcf_file), out_pos_name, '-p', 'isec', '-O', 'z']
    context.runner.call(job, cmd, work_dir=work_dir)

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', 'isec/0000.vcf.gz'], work_dir=work_dir)

    neg_control_vcf_id = context.write_output_file(job, os.path.join(work_dir, 'isec', '0000.vcf.gz'),
                                                   out_store_path = out_neg_name)
    neg_control_tbi_id = context.write_output_file(job, os.path.join(work_dir, 'isec', '0000.vcf.gz.tbi'),
                                                   out_store_path = out_neg_name + '.tbi')

    return pos_control_vcf_id, pos_control_tbi_id, neg_control_vcf_id, neg_control_tbi_id
    
    
def construct_main(context, options):
    """
    Wrapper for vg constructing. 
    """
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            inputFastaFileID = toil.importFile(options.fasta)
            if options.vcf:
                inputVCFFileID = toil.importFile(options.vcf)
                inputVCFName = os.path.basename(options.vcf)
                inputTBIFileID = toil.importFile(options.vcf + '.tbi')
            else:
                inputVCFFileID = None
                inputVCFName = None
                inputTBIFileID = None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_construct_with_controls, context, options.control_sample,
                                     inputFastaFileID,
                                     os.path.basename(options.fasta),
                                     inputVCFFileID, inputVCFName,
                                     inputTBIFileID,
                                     options.max_node_size,
                                     options.alt_paths,
                                     options.flat_alts,
                                     options.regions,
                                     merge_output_name = options.graph_name)
            # Run the job
            toil.start(root_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
