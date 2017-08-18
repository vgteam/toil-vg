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
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_index import run_xg_indexing, run_gcsa_prep
logger = logging.getLogger(__name__)

# from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/Illumina_PlatinumGenomes_NA12877_NA12878_09162015/IlluminaPlatinumGenomes-user-guide.pdf
CEPH_SAMPLES="NA12889 NA12890 NA12891 NA12892 NA12877 NA12878 NA12879 NA12880 NA12881 NA12882 NA12883 NA12884 NA12885 NA12886 NA12887 NA12888 NA12893".split()

def construct_subparser(parser):
    """
    Create a subparser for construction.  Should pass in results of subparsers.add_parser()
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
    parser.add_argument("--max_node_size", type=int, default=32,
                        help="Maximum node length")
    parser.add_argument("--alt_paths", action="store_true",
                        help="Save paths for alts with variant ID (required for GPBWT)")
    parser.add_argument("--flat_alts", action="store_true",
                        help="flat alts")
    parser.add_argument("--control_sample", default=None,
                        help="Make a positive and negative control using this sample")
    parser.add_argument("--construct_cores", type=int,
                        help="Number of cores for vg construct")
    parser.add_argument("--out_name", required=True,
                        help="Name used for output graphs and indexes")
    parser.add_argument("--merge_graphs", action="store_true",
                        help="Merge all regions into one graph")
    parser.add_argument("--filter_ceph", action="store_true",
                        help="Filter out all variants specific to the CEPH pedigree, which includes NA12878")
    parser.add_argument("--filter_samples", nargs='+',
                        help="Filter out all variants specific to given samples")
    parser.add_argument("--gcsa_index", action="store_true",
                        help="Make a gcsa index for each output graph")
    parser.add_argument("--xg_index", action="store_true",
                        help="Make an xg index for each output graph")
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def run_generate_input_vcfs(job, context, sample, fasta_id, vcf_id, vcf_name, tbi_id,
                            regions, output_name, filter_samples = None):
    """
    Preprocessing step to make a bunch of vcfs if wanted:
    - positive control
    - negative control
    - family filter
    returns a list of (vcf_id, vcf_name, tbi_id, merge_name, region_names) tuples
    """
    # our input vcf 
    output = [(vcf_id, vcf_name, tbi_id, output_name,
               [c.replace(':','-') for c in regions] if regions else None)]
    
    # our positive and negative controls
    if sample:
        make_controls = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id, sample,
                                          cores=context.config.construct_cores,
                                          memory=context.config.construct_mem,
                                          disk=context.config.construct_disk)
        pos_control_vcf_id, pos_control_tbi_id = make_controls.rv(0), make_controls.rv(1)
        neg_control_vcf_id, neg_control_tbi_id = make_controls.rv(2), make_controls.rv(3)

        vcf_base = os.path.basename(vcf_name.rstrip('.gz').rstrip('.vcf'))
        pos_control_vcf_name = '{}_{}.vcf.gz'.format(vcf_base, sample)
        neg_control_vcf_name = '{}_minus_{}.vcf.gz'.format(vcf_base, sample)
        if regions:
            pos_region_names = [output_name + '_' + c.replace(':','-') + '_{}'.format(sample) for c in regions]
            neg_region_names = [output_name + '_' + c.replace(':','-') + '_minus_{}'.format(sample) for c in regions]
        else:
            pos_region_names = None
            neg_region_names = None
        pos_output_name = output_name.rstrip('.vg') + '_{}.vg'.format(sample)
        neg_output_name = output_name.rstrip('.vg') + '_minus_{}.vg'.format(sample)

        output.append((pos_control_vcf_id, pos_control_vcf_name, pos_control_tbi_id,
                       pos_output_name, pos_region_names))

        output.append((neg_control_vcf_id, neg_control_vcf_name, neg_control_tbi_id,
                       neg_output_name, neg_region_names))

    # our family filter
    if filter_samples:
        filter_job = job.addChildJobFn(run_filter_vcf_samples, context, vcf_id, vcf_name, tbi_id,
                                       filter_samples,
                                       cores=context.config.construct_cores,
                                       memory=context.config.construct_mem,
                                       disk=context.config.construct_disk)

        filter_vcf_id, filter_tbi_id = filter_job.rv(0), filter_job.rv(1)

        vcf_base = os.path.basename(vcf_name.rstrip('.gz').rstrip('.vcf'))
        filter_vcf_name = '{}_filter.vcf.gz'.format(vcf_base)
        if regions:
            filter_region_names = [output_name + '_' + c.replace(':','-') + '_filter' for c in regions]
        else:
            filter_region_names = None
        filter_output_name = output_name.rstrip('.vg') + '_filter.vg'

        output.append((filter_vcf_id, filter_vcf_name, filter_tbi_id,
                       filter_output_name, filter_region_names))

    return output

    
def run_construct_all(job, context, fasta_id, fasta_name, vcf_inputs, 
                      max_node_size, alt_paths, flat_alts, regions, merge_graphs,
                      sort_ids, join_ids, gcsa_index, xg_index):
    """ 
    construct many graphs in parallel, optionally doing indexing too.  vcf_inputs
    is a list of tuples as created by  run_generate_input_vcfs
    """

    output = []
    
    for (vcf_id, vcf_name, tbi_id, output_name, region_names) in vcf_inputs:
        merge_output_name = output_name if merge_graphs or not regions or len(regions) < 2 else None
        construct_job = job.addChildJobFn(run_construct_genome_graph, context, fasta_id,
                                          fasta_name, vcf_id, vcf_name, tbi_id,
                                          max_node_size, alt_paths, flat_alts, regions,
                                          region_names, sort_ids, join_ids, merge_output_name)
        vg_ids = construct_job.rv()
        vg_names = output_name if merge_graphs else region_names

        if gcsa_index:
            gcsa_job = job.addFollowOnJobFn(run_gcsa_prep, context, vg_ids,
                                            vg_names, output_name, regions)
            gcsa_id = gcsa_job.rv(0)
            lcp_id = gcsa_job.rv(1)
        else:
            gcsa_id = None
            lcp_id = None
            
        if xg_index:
            xg_id = job.addFollowOnJobFn(run_xg_indexing, context, vg_ids,
                                         vg_names, output_name).rv()
        else:
            xg_id = None

        output.append((vg_ids, vg_names, gcsa_id, lcp_id, xg_id))
    return output
                

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
                                                  sort_ids=sort_ids,
                                                  cores=context.config.construct_cores,
                                                  memory=context.config.construct_mem,
                                                  disk=context.config.construct_disk).rv())

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

    if merge_output_name and not merge_output_name.endswith('.vg'):
        merge_output_name += '.vg'

    # if there's nothing to do, just write the files and return
    if len(region_graph_ids) == 1 or not (join_ids or merge_output_name):
        out_ids = []
        for region_file in region_files:
            out_ids.append(context.write_output_file(job, os.path.join(work_dir, region_file),
                                                     out_store_path = merge_output_name))
        return out_ids

    if join_ids:
        # join the ids
        cmd = ['vg', 'ids', '-j'] + region_files
        context.runner.call(job, cmd, work_dir=work_dir)

    if merge_output_name is not None:
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

    cmd = ['vg', 'construct', '-r', os.path.basename(fasta_file)]
    if vcf_id:
        cmd += ['-v', os.path.basename(vcf_file)]
    if region:
        cmd += ['-R', region]
    if is_chrom:
        cmd += ['-C']
    if max_node_size:
        cmd += ['-m', max_node_size]
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

def run_filter_vcf_samples(job, context, vcf_id, vcf_name, tbi_id, samples):
    """ Use vcflib to remove all variants specifc to a set of samples.
    
    This is extremely slow.  Will want to parallelize if doing often on large VCFs
    (or rewrite custom tool?  I think running time sunk in vcffixup recomputing allele freqs
    which is overkill)
    """
    if not samples:
        return vcf_id, tbi_id
    
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    cmd = [['vcfremovesamples', os.path.basename(vcf_file)] + samples]
    cmd.append(['vcffixup', '-'])
    cmd.append(['vcffilter', '-f', 'AC > 0'])

    vcf_base = os.path.basename(vcf_name.rstrip('.gz').rstrip('.vcf'))
    filter_vcf_name = '{}_filter.vcf'.format(vcf_base)

    with open(os.path.join(work_dir, filter_vcf_name), 'w') as out_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile=out_file)

    # bgzip in separate command because docker interface requires (big waste of time/space)
    # note: tried to use Bio.bgzf.open above to get around but it doesn't seem to work
    # with streaming
    context.runner.call(job, ['bgzip', filter_vcf_name], work_dir=work_dir)
    filter_vcf_name += '.gz'

    out_vcf_id = context.write_output_file(job, os.path.join(work_dir, filter_vcf_name))

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', filter_vcf_name],
                        work_dir=work_dir)
                                        
    out_tbi_id = context.write_output_file(job, os.path.join(work_dir, filter_vcf_name) + '.tbi')
    
    return out_vcf_id, out_tbi_id
    
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
    
    # remove anything that's not alt (probably cleaner way to do this)
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

    # Merge up all filter samples
    filter_samples = []
    if options.filter_samples:
        filter_samples += options.filter_samples
    if options.filter_ceph:
        filter_samples += CEPH_SAMPLES
    filter_samples = list(set(filter_samples))

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

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context)

            # Automatically make and name a bunch of vcfs
            vcf_job = init_job.addFollowOnJobFn(run_generate_input_vcfs, context, options.control_sample,
                                                inputFastaFileID, 
                                                inputVCFFileID, inputVCFName,
                                                inputTBIFileID, 
                                                options.regions,
                                                options.out_name,
                                                filter_samples)

            # Cosntruct graphs
            vcf_job.addFollowOnJobFn(run_construct_all, context, inputFastaFileID,
                                     os.path.basename(options.fasta), vcf_job.rv(),
                                     options.max_node_size, options.alt_paths,
                                     options.flat_alts, options.regions, options.merge_graphs,
                                     True, True, options.gcsa_index, options.xg_index)
            
            # Run the workflow
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
