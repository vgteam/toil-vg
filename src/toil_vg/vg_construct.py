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
import os

from math import ceil
from subprocess import Popen, PIPE

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import *
from toil_vg.context import Context, run_write_info_to_outstore
from toil_vg.vg_index import run_xg_indexing, run_indexing, index_parse_args, index_toggle_parse_args
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

    parser.add_argument("--fasta", required=True, type=make_url, nargs='+',
                        help="Reference sequence in fasta or fasta.gz (single fasta or 1/region in same order as --regions)")
    parser.add_argument("--vcf", default=None, type=make_url, nargs='+',
                        help="Variants to make graph from (single vcf or 1/region in same order as --regions)")
    parser.add_argument("--regions", nargs='+',
                        help="1-based inclusive VCF coordinates in the form of SEQ or SEQ:START-END")
    parser.add_argument("--fasta_regions", action="store_true",
                        help="Infer regions from fasta file.  If multiple vcfs specified, any regions found that are not in --regions will be added without variants (useful for decoy sequences)")    
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
    parser.add_argument("--haplo_sample", type=str,
                        help="Make haplotype thread graphs (for simulating from) for this sample")
    parser.add_argument("--primary", action="store_true",
                        help="Make the primary graph (no variants)")
    parser.add_argument("--no_base", action="store_true",
                        help="Do not construct base graph from input vcf.  Only make optional controls")
    parser.add_argument("--min_af", type=float, default=None,
                        help="Create a control using the given minium allele frequency")

    # Add common indexing options shared with vg_index
    index_toggle_parse_args(parser)
    index_parse_args(parser)

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def validate_construct_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(not options.haplo_sample or options.regions,
            '--regions required with --haplo_sample')
    require(not (options.haplo_sample and options.merge_graphs),
            '--merge_graphs not currently supported with --haplo_sample')
    require(not (options.haplo_sample and options.fasta_regions),
            '--fasta_regions not currently supported with --haplo_sample')    
    require(not options.vcf or len(options.vcf) == 1 or not options.regions or
            len(options.vcf) <= len(options.regions),
            'if many vcfs specified, cannot have more vcfs than --regions')
    require(len(options.fasta) == 1 or len(options.fasta) == len(options.regions),
            'if many fastas specified, must be same number as --regions')
    require(len(options.fasta) == 1 or not options.fasta_regions,
            '--fasta_regions currently only works when single fasta specified with --fasta')
    require(not options.gbwt_index or options.xg_index,
            '--xg_index required with --gbwt_index')
    # TODO: It seems like somne of this code is designed to run multiple regions
    # in parallel, but the indexing code always indexes them together.
    require(not options.gbwt_index or len(options.vcf) >= 1,
            '--gbwt_index requires --vcf')
            
    
def run_unzip_fasta(job, context, fasta_id, fasta_name):
    """
    vg construct doesn't work with zipped fasta, so we run this on input fasta that end in .gz
    """
    
    work_dir = job.fileStore.getLocalTempDir()

    # Download input files
    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    job.fileStore.readGlobalFile(fasta_id, fasta_file, mutable=True)
    context.runner.call(job, ['bgzip', '-d', os.path.basename(fasta_file)], work_dir=work_dir)

    return context.write_intermediate_file(job, fasta_file[:-3])

def run_scan_fasta_sequence_names(job, context, fasta_id, fasta_name, regions = None):
    """
    if no regions specified, scrape them out of the (uncompressed) fasta
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Download input files
    fasta_file = os.path.join(work_dir, os.path.basename(fasta_name))
    job.fileStore.readGlobalFile(fasta_id, fasta_file)
    
    # reluctant to use slow python library, so just running grep instead
    cmd = ['grep', '>', os.path.basename(fasta_file)]
    grep_output = context.runner.call(job, cmd, work_dir = work_dir,
                                      check_output = True, tool_name='bgzip')

    # just taking first whitespace-separated token.  that's what corresponds to hs37d5 vcf
    seq_names = [] if not regions else regions
    for line in grep_output.split('\n'):
        if len(line) > 1:
            name = line.split()[0]
            if name.startswith('>') and (not regions or name[1:] not in regions):
                seq_names.append(name[1:])
    
    return seq_names    
        
def run_generate_input_vcfs(job, context, sample, vcf_ids, vcf_names, tbi_ids,
                            regions, output_name, filter_samples = None,
                            haplo_sample = None, do_primary = False, min_af = None,
                            make_base_graph = True):
    """
    Preprocessing step to make a bunch of vcfs if wanted:
    - positive control
    - negative control
    - family filter
    - primary
    - thresholded by a given minimum allele frequency
    returns a dictionary of name -> (vcf_id, vcf_name, tbi_id, merge_name, region_names) tuples
    where name can be used to, ex, tell the controls apart
    """
    # our input vcf
    output = dict()
    if make_base_graph:
        output[output_name] = [vcf_ids, vcf_names, tbi_ids, output_name,
                               [output_name + '_' + c.replace(':','-') for c in regions] if regions else None]
    
    # our positive and negative controls
    if sample:
        pos_control_vcf_ids, pos_control_tbi_ids = [], []
        neg_control_vcf_ids, neg_control_tbi_ids = [], []
        pos_control_vcf_names, neg_control_vcf_names = [], []
        
        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            make_controls = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id, sample,
                                              cores=context.config.construct_cores,
                                              memory=context.config.construct_mem,
                                              disk=context.config.construct_disk)
            pos_control_vcf_ids.append(make_controls.rv(0))
            pos_control_tbi_ids.append(make_controls.rv(1))
            neg_control_vcf_ids.append(make_controls.rv(2))
            neg_control_tbi_ids.append(make_controls.rv(3))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            pos_control_vcf_names.append('{}_{}.vcf.gz'.format(vcf_base, sample))
            neg_control_vcf_names.append('{}_minus_{}.vcf.gz'.format(vcf_base, sample))
        if regions:
            pos_region_names = [output_name + '_{}'.format(sample)  + '_' + c.replace(':','-') for c in regions]
            neg_region_names = [output_name + '_minus_{}'.format(sample) + '_' + c.replace(':','-')  for c in regions]
        else:
            pos_region_names = None
            neg_region_names = None
        pos_output_name = remove_ext(output_name, '.vg') + '_{}.vg'.format(sample)
        neg_output_name = remove_ext(output_name, '.vg') + '_minus_{}.vg'.format(sample)

        output['pos-control'] = [pos_control_vcf_ids, pos_control_vcf_names, pos_control_tbi_ids,
                                 pos_output_name, pos_region_names]

        output['neg-control'] = [neg_control_vcf_ids, neg_control_vcf_names, neg_control_tbi_ids,
                                 neg_output_name, neg_region_names]

        if haplo_sample and haplo_sample == sample:
            output['haplo'] = output['pos-control']

    # our family filter
    if filter_samples:
        filter_vcf_ids, filter_tbi_ids = [], []
        filter_vcf_names = []

        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            filter_job = job.addChildJobFn(run_filter_vcf_samples, context, vcf_id, vcf_name, tbi_id,
                                           filter_samples,
                                           cores=context.config.construct_cores,
                                           memory=context.config.construct_mem,
                                           disk=context.config.construct_disk)

            filter_vcf_ids.append(filter_job.rv(0))
            filter_tbi_ids.append(filter_job.rv(1))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            filter_vcf_names.append('{}_filter.vcf.gz'.format(vcf_base))
        if regions:
            filter_region_names = [output_name + '_filter'  + '_' + c.replace(':','-') for c in regions]
        else:
            filter_region_names = None
        filter_output_name = remove_ext(output_name, '.vg') + '_filter.vg'

        output['filter'] = [filter_vcf_ids, filter_vcf_names, filter_tbi_ids,
                            filter_output_name, filter_region_names]

    # we want a vcf to make a gpbwt out of for making haplo graphs
    # we re-use the vcf from the positive control if available, but we give it a
    # different name and construct different .vg graphs going forward to allow for
    # different construction (ie the haplo graph will always get alt paths that we don't
    # necessarily want in the control)
    if haplo_sample:
        hap_control_vcf_ids, hap_control_tbi_ids = [], []
        hap_control_vcf_names = []

        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            if haplo_sample != sample:
                make_controls = job.addChildJobFn(run_make_control_vcfs, context, vcf_id, vcf_name, tbi_id, haplo_sample,
                                                  pos_only = True,
                                                  cores=context.config.construct_cores,
                                                  memory=context.config.construct_mem,
                                                  disk=context.config.construct_disk)
                hap_control_vcf_ids.append(make_controls.rv(0))
                hap_control_tbi_ids.append(make_controls.rv(1))
            else:
                hap_control_vcf_ids = pos_control_vcf_ids
                hap_control_tbi_ids = pos_control_tbi_ids
            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            hap_control_vcf_names.append('{}_{}_haplo.vcf.gz'.format(vcf_base, haplo_sample))
        if regions:
            hap_region_names = [output_name + '_{}'.format(haplo_sample)  + '_' + c.replace(':','-') for c in regions]
        else:
            hap_region_names = None
        hap_output_name = remove_ext(output_name, '.vg') + '_{}_haplo.vg'.format(haplo_sample)
        
        output['haplo'] = [hap_control_vcf_ids, hap_control_vcf_names, hap_control_tbi_ids,
                           hap_output_name, hap_region_names]

    if do_primary:
        if regions:
            primary_region_names = ['primary'  + '_' + c.replace(':','-') for c in regions]
        else:
            primary_region_names = None

        primary_output_name = 'primary.vg' if '_' not in output_name else 'primary' + output_name[output_name.find('_')+1:]
        output['primary'] = [[], [], [], primary_output_name, primary_region_names]

    if min_af is not None:
        af_vcf_ids, af_tbi_ids = [], []
        af_vcf_names = []

        for vcf_id, vcf_name, tbi_id in zip(vcf_ids, vcf_names, tbi_ids):
            af_job = job.addChildJobFn(run_min_allele_filter_vcf_samples, context, vcf_id, vcf_name, tbi_id,
                                       min_af,
                                       cores=context.config.construct_cores,
                                       memory=context.config.construct_mem,
                                       disk=context.config.construct_disk)

            af_vcf_ids.append(af_job.rv(0))
            af_tbi_ids.append(af_job.rv(1))

            vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
            af_vcf_names.append('{}_minaf_{}.vcf.gz'.format(vcf_base, min_af))
        if regions:
            af_region_names = [output_name + '_af' + '_' + c.replace(':','-') for c in regions]
        else:
            af_region_names = None
        af_output_name = remove_ext(output_name, '.vg') + '_minaf_{}.vg'.format(min_af)

        output['minaf'] = [af_vcf_ids, af_vcf_names, af_tbi_ids,
                           af_output_name, af_region_names]

    # pad out vcf lists with nones so they are the same size as regions
    # since we allow fasta regions that dont have corresponding vcf
    # note: single vcf , multiple regions case not handled here as it's
    # treated below (the same vcf is given to each region)
    if regions and len(regions) > len(vcf_ids) and len(vcf_ids) != 1:
        padding = [None] * (len(regions) - len(vcf_ids))
        for key, val in output.items():
            val[0] += padding
            val[1] += padding
            val[2] += padding

    return output

    
def run_construct_all(job, context, fasta_ids, fasta_names, vcf_inputs, 
                      max_node_size, alt_paths, flat_alts, regions,
                      merge_graphs = False, sort_ids = False, join_ids = False,
                      gcsa_index = False, xg_index = False, gbwt_index = False, snarls_index = False,
                      haplo_sample = None, haplotypes = [0,1], gbwt_prune = False):
    """ 
    construct many graphs in parallel, optionally doing indexing too. vcf_inputs
    is a list of tuples as created by run_generate_input_vcfs
    
    Returns a list of tuples of the form (vg_ids, vg_names, indexes), where
    indexes is the index dict from index type to file ID.
    
    """

    output = []
    
    for name, (vcf_ids, vcf_names, tbi_ids, output_name, region_names) in vcf_inputs.items():
        merge_output_name = output_name if merge_graphs or not regions or len(regions) < 2 else None
        output_name_base = remove_ext(output_name, '.vg')
        gpbwt = name == 'haplo'
        construct_job = job.addChildJobFn(run_construct_genome_graph, context, fasta_ids,
                                          fasta_names, vcf_ids, vcf_names, tbi_ids,
                                          max_node_size, gbwt_index or gpbwt or alt_paths,
                                          flat_alts, regions,
                                          region_names, sort_ids, join_ids, name, merge_output_name)

        vg_ids = construct_job.rv(0)
        mapping_id = construct_job.rv(1)
        vg_names = [merge_output_name] if merge_graphs or not regions or len(regions) < 2 else region_names

        vg_names = [remove_ext(i, '.vg') + '.vg' for i in vg_names]

        if not regions:
            chroms = []
        else:
            chroms = [p.split(':')[0] for p in regions]

        # strip nones out of vcf list            
        input_vcf_ids = []
        input_tbi_ids = []
        if gpbwt or gbwt_index:
            for vcf_id, tbi_id in zip(vcf_ids, tbi_ids):
                if vcf_id and tbi_id:
                    input_vcf_ids.append(vcf_id)
                    input_tbi_ids.append(tbi_id)
                else:
                    assert vcf_id == None and tbi_id == None
        
        indexing_job = construct_job.addFollowOnJobFn(run_indexing, context, vg_ids,
                                                      vg_names, output_name_base, chroms,
                                                      input_vcf_ids, input_tbi_ids,
                                                      node_mapping_id = mapping_id,
                                                      skip_xg=not xg_index, skip_gcsa=not gcsa_index,
                                                      skip_id_ranges=True, skip_snarls=not snarls_index,
                                                      make_gbwt=gbwt_index, gbwt_prune=gbwt_prune)
        indexes = indexing_job.rv()    

        if gpbwt:
            assert(haplo_sample is not None)
            haplo_job = construct_job.addFollowOnJobFn(run_make_haplo_graphs, context,
                                                       input_vcf_ids, input_tbi_ids,
                                                       vcf_names, vg_ids, vg_names, output_name_base, regions,
                                                       haplo_sample, haplotypes)

            # we want an xg index from our thread graphs to pass to vg sim for each haplotype
            for haplotype in haplotypes:
                haplo_xg_job = haplo_job.addFollowOnJobFn(run_xg_indexing, context, haplo_job.rv(haplotype),
                                                          vg_names,
                                                          output_name_base + '_thread_{}'.format(haplotype),
                                                          cores=context.config.xg_index_cores,
                                                          memory=context.config.xg_index_mem,
                                                          disk=context.config.xg_index_disk)
                                                          
                # TODO: file IDs are not exposed to caller
    
    

        output.append((vg_ids, vg_names, indexes))
    return output
                

def run_construct_genome_graph(job, context, fasta_ids, fasta_names, vcf_ids, vcf_names, tbi_ids,
                               max_node_size, alt_paths, flat_alts, regions, region_names,
                               sort_ids, join_ids, name, merge_output_name):
    """ construct graph(s) from several regions in parallel.  we could eventually generalize this
    to accept multiple vcfs and/or fastas if needed, as well as to determine regions from file,
    but for now we only accept single files, and require region list.
    """

    # encapsulate follow-on
    child_job = Job()
    job.addChild(child_job)
    
    work_dir = job.fileStore.getLocalTempDir()

    if not regions:
        regions, region_names = [None], ['genome']        

    region_graph_ids = []    
    for i, (region, region_name) in enumerate(zip(regions, region_names)):
        vcf_id = None if not vcf_ids or i >= len(vcf_ids) else vcf_ids[0] if len(vcf_ids) == 1 else vcf_ids[i]
        vcf_name = None if not vcf_names or i >= len(vcf_names) else vcf_names[0] if len(vcf_names) == 1 else vcf_names[i]
        tbi_id = None if not tbi_ids or i >= len(tbi_ids) else tbi_ids[0] if len(tbi_ids) == 1 else tbi_ids[i]
        fasta_id = fasta_ids[0] if len(fasta_ids) == 1 else fasta_ids[i]
        fasta_name = fasta_names[0] if len(fasta_names) == 1 else fasta_names[i]
        region_graph_ids.append(child_job.addChildJobFn(run_construct_region_graph, context,
                                                  fasta_id, fasta_name,
                                                  vcf_id, vcf_name, tbi_id, region, region_name,
                                                  max_node_size, alt_paths, flat_alts,
                                                  # todo: bump as command line option?
                                                  #       also, needed if we update vg docker image?
                                                  is_chrom=not region or ':' not in region,
                                                  sort_ids=sort_ids,
                                                  cores=context.config.construct_cores,
                                                  memory=context.config.construct_mem,
                                                  disk=context.config.construct_disk).rv())

    return child_job.addFollowOnJobFn(run_join_graphs, context, region_graph_ids, join_ids,
                                      region_names, name, merge_output_name).rv()


def run_join_graphs(job, context, region_graph_ids, join_ids, region_names, name, merge_output_name = None):
    """
    join the ids of some graphs.  if a merge_output_name is given, cat them all together as well
    this function saves output to the outstore.  also does the node mapping file 
    """
        
    work_dir = job.fileStore.getLocalTempDir()

    # download graph for each region
    region_files = []
    for region_graph_id, region_name in zip(region_graph_ids, region_names):
        region_file = '{}.vg'.format(region_name)
        job.fileStore.readGlobalFile(region_graph_id, os.path.join(work_dir, region_file), mutable=True)
        region_files.append(region_file)

    if merge_output_name:
        merge_output_name = remove_ext(merge_output_name, '.vg') + '.vg'

    # if there's nothing to do, just write the files and return
    if not join_ids:
        out_ids = []
        for region_file in region_files:
            out_ids.append(context.write_output_file(job, os.path.join(work_dir, region_file),
                                                     out_store_path = merge_output_name))
        return out_ids, None

    mapping_file = merge_output_name[:-3] if merge_output_name else name
    mapping_file = os.path.join(work_dir, mapping_file + '.mapping')

    if join_ids:
        # join the ids
        cmd = ['vg', 'ids', '--join', '--mapping', os.path.basename(mapping_file)] + region_files
        context.runner.call(job, cmd, work_dir=work_dir)

    if merge_output_name is not None:
        assert merge_output_name[:-3] not in region_names
        with open(os.path.join(work_dir, merge_output_name), 'w') as merge_file:
            for region_file in region_files:
                with open(os.path.join(work_dir, region_file)) as cf:
                    shutil.copyfileobj(cf, merge_file)
        out_graphs = [context.write_output_file(job, os.path.join(work_dir, merge_output_name))]
    else:
        out_graphs = [context.write_output_file(job, os.path.join(work_dir, f)) for f in region_files]

    mapping_id = context.write_output_file(job, mapping_file)
    return out_graphs, mapping_id
        
    
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

    cmd = ['vg', 'construct', '--reference', os.path.basename(fasta_file)]
    if vcf_id:
        cmd += ['--vcf', os.path.basename(vcf_file)]
    if region:
        cmd += ['--region', region]
        if is_chrom:
            cmd += ['--region-is-chrom']
    if max_node_size:
        cmd += ['--node-max', max_node_size]
    if alt_paths:
        cmd += ['--alt-paths']
    if job.cores:
        cmd += ['--threads', job.cores]

    if sort_ids:
        cmd = [cmd, ['vg', 'ids', '--sort', '-']]

    vg_path = os.path.join(work_dir, region_name)
    with open(vg_path, 'w') as vg_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile = vg_file)

    return context.write_intermediate_file(job, vg_path)

def run_filter_vcf_samples(job, context, vcf_id, vcf_name, tbi_id, samples):
    """ Use vcflib to remove all variants specifc to a set of samples.
    """
    if not samples:
        return vcf_id, tbi_id
    
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    # Warning: This VCF is only going to have the samples we want to filter out.
    # Will be okay for constructing graphs, but cause problems for anything that
    # needs sample information. 
    cmd = ['bcftools', 'view', os.path.basename(vcf_file), '--exclude-private',
           '--samples', ','.join(samples), '--force-samples', '--output-type', 'z']

    vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
    filter_vcf_name = '{}_filter.vcf.gz'.format(vcf_base)

    with open(os.path.join(work_dir, filter_vcf_name), 'w') as out_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile=out_file)

    out_vcf_id = context.write_output_file(job, os.path.join(work_dir, filter_vcf_name))

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', filter_vcf_name],
                        work_dir=work_dir)
                                        
    out_tbi_id = context.write_output_file(job, os.path.join(work_dir, filter_vcf_name) + '.tbi')
    
    return out_vcf_id, out_tbi_id
    
def run_make_control_vcfs(job, context, vcf_id, vcf_name, tbi_id, sample, pos_only = False):
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

    out_pos_name = remove_ext(remove_ext(os.path.basename(vcf_name), '.gz'), '.vcf')
    out_neg_name = out_pos_name + '_minus_{}.vcf.gz'.format(sample)
    out_pos_name += '_{}.vcf.gz'.format(sample)

    with open(os.path.join(work_dir, out_pos_name), 'w') as out_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile = out_file)

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', out_pos_name], work_dir=work_dir)

    pos_control_vcf_id = context.write_output_file(job, os.path.join(work_dir, out_pos_name))
    pos_control_tbi_id = context.write_output_file(job, os.path.join(work_dir, out_pos_name + '.tbi'))

    if pos_only:
        return pos_control_vcf_id, pos_control_tbi_id, None, None

    # subtract the positive control to make the negative control
    cmd = ['bcftools', 'isec', os.path.basename(vcf_file), out_pos_name, '-p', 'isec', '-O', 'z']
    context.runner.call(job, cmd, work_dir=work_dir)

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', 'isec/0000.vcf.gz'], work_dir=work_dir)

    neg_control_vcf_id = context.write_output_file(job, os.path.join(work_dir, 'isec', '0000.vcf.gz'),
                                                   out_store_path = out_neg_name)
    neg_control_tbi_id = context.write_output_file(job, os.path.join(work_dir, 'isec', '0000.vcf.gz.tbi'),
                                                   out_store_path = out_neg_name + '.tbi')

    return pos_control_vcf_id, pos_control_tbi_id, neg_control_vcf_id, neg_control_tbi_id

def run_min_allele_filter_vcf_samples(job, context, vcf_id, vcf_name, tbi_id, min_af):
    """
    filter a vcf by allele frequency using bcftools --min-af
    """
    if not min_af:
        return vcf_id, tbi_id
    
    work_dir = job.fileStore.getLocalTempDir()

    vcf_file = os.path.join(work_dir, os.path.basename(vcf_name))
    job.fileStore.readGlobalFile(vcf_id, vcf_file)
    job.fileStore.readGlobalFile(tbi_id, vcf_file + '.tbi')

    vcf_base = os.path.basename(remove_ext(remove_ext(vcf_name, '.gz'), '.vcf'))
    af_vcf_name = '{}_minaf_{}.vcf.gz'.format(vcf_base, min_af)

    cmd = ['bcftools', 'view', '--min-af', min_af, '-O', 'z', os.path.basename(vcf_file)]
    with open(os.path.join(work_dir, af_vcf_name), 'w') as out_file:
        context.runner.call(job, cmd, work_dir = work_dir, outfile=out_file)

    out_vcf_id = context.write_output_file(job, os.path.join(work_dir, af_vcf_name))

    context.runner.call(job, ['tabix', '-f', '-p', 'vcf', af_vcf_name],
                        work_dir=work_dir)
                                        
    out_tbi_id = context.write_output_file(job, os.path.join(work_dir, af_vcf_name) + '.tbi')
    
    return out_vcf_id, out_tbi_id

def run_make_haplo_graphs(job, context, vcf_ids, tbi_ids, vcf_names, vg_ids, vg_names,
                          output_name, regions, sample, haplotypes):
    """
    make some haplotype graphs for threads in a gpbwt.  regions must be defined since we use the
    chromosome name to get the threads
    """

    assert(sample is not None)

    # ith element will be a list of graphs (1 list / region) for haplotype i
    thread_vg_ids = []
    for h in haplotypes:
        thread_vg_ids.append([])

    # make sure we're only dealing with chrom names (should probably be error otherwise)
    chroms = [region[0:region.find(':')] if ':' in region else region for region in regions]

    # validate options should enforce this but check to be sure assumptions met to avoid
    # returning nonsense
    assert len(vg_ids) == len(regions)
    assert len(vcf_ids) == 1 or len(vcf_ids) == len(regions)
    assert len(tbi_ids) == len(vcf_ids)
    
    logger.info('Making haplo graphs for chromosomes {}'.format(chroms))
    
    for i, (vg_id, vg_name, region) in enumerate(zip(vg_ids, vg_names, chroms)):
        vcf_name = vcf_names[0] if len(vcf_names) == 1 else vcf_names[i]
        vcf_id = vcf_ids[0] if len(vcf_names) == 1 else vcf_ids[i]
        tbi_id = tbi_ids[0] if len(vcf_names) == 1 else tbi_ids[i]
            
        # index the graph and vcf to make the gpbwt
        gpbwt_name = '{}-gpbwt'.format(remove_ext(vg_name, '.vg'))
        gpbwt_job = job.addChildJobFn(run_xg_indexing, context, [vg_id], [vg_name],
                                      gpbwt_name, vcf_id, tbi_id,
                                      cores=context.config.xg_index_cores,
                                      memory=context.config.xg_index_mem,
                                      disk=context.config.xg_index_disk)
        xg_id = gpbwt_job.rv(0)
        

        # make a thread graph from the xg
        hap_job = gpbwt_job.addFollowOnJobFn(run_make_haplo_thread_graphs, context, vg_id, vg_name,
                                             output_name, [region], xg_id, sample, haplotypes,
                                             cores=context.config.construct_cores,
                                             memory=context.config.construct_mem,
                                             disk=context.config.construct_disk)
        for j in range(len(haplotypes)):
            thread_vg_ids[j].append(hap_job.rv(j))

    return thread_vg_ids

def run_make_haplo_thread_graphs(job, context, vg_id, vg_name, output_name, chroms, xg_id,
                                sample, haplotypes):
    """
    make some haplotype graphs for threads in a gpbwt
    """
    work_dir = job.fileStore.getLocalTempDir()

    xg_path = os.path.join(work_dir, vg_name[:-3] + '-gpbwt.xg')
    job.fileStore.readGlobalFile(xg_id, xg_path)

    vg_path = os.path.join(work_dir, vg_name)
    job.fileStore.readGlobalFile(vg_id, vg_path)

    thread_vg_ids = []

    for hap in haplotypes:
        
        # This can't work if the sample is None and we want any haplotypes
        assert(sample is not None)

        try:

            tag = '_{}'.format(chroms[0]) if len(chroms) == 1 else ''
            thread_path = os.path.join(work_dir, '{}{}_thread_{}_merge.vg'.format(output_name, tag, hap))
            logger.info('Creating thread graph {}'.format(thread_path))
            with open(thread_path, 'w') as thread_file:
                # strip paths from our original graph            
                cmd = ['vg', 'mod', '-D', os.path.basename(vg_path)]
                context.runner.call(job, cmd, work_dir = work_dir, outfile = thread_file)

                # get haplotype thread paths from the gpbwt
                cmd = ['vg', 'find', '-x', os.path.basename(xg_path)]
                for chrom in chroms:
                    cmd += ['-q', '_thread_{}_{}_{}'.format(sample, chrom, hap)]
                context.runner.call(job, cmd, work_dir = work_dir, outfile = thread_file)
                
            thread_path_trim = os.path.join(work_dir, '{}{}_thread_{}.vg'.format(output_name, tag, hap))
            logger.info('Creating trimmed thread graph {}'.format(thread_path_trim))
            with open(thread_path_trim, 'w') as thread_file:
                # Then we trim out anything other than our thread path
                cmd = [['vg', 'mod', '-N', os.path.basename(thread_path)]]
                # And get rid of our thread paths since they take up lots of space when re-indexing
                filter_cmd = ['vg', 'mod', '-']
                for chrom in chroms:
                    filter_cmd += ['-r', chrom]
                cmd.append(filter_cmd)
                context.runner.call(job, cmd, work_dir = work_dir, outfile = thread_file)

            thread_vg_ids.append(context.write_output_file(job, thread_path_trim))
            
        except:
            # Dump everything we need to replicate the thread extraction
            logging.error("Thread extraction failed. Dumping files.")

            context.write_output_file(job, vg_path)
            context.write_output_file(job, xg_path)
            
            raise

    logger.info('Got {} thread file IDs'.format(len(thread_vg_ids)))

    return thread_vg_ids

def construct_main(context, options):
    """
    Wrapper for vg constructing. 
    """

    # check some options
    validate_construct_options(options)

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
            inputFastaFileIDs = [toil.importFile(fasta) for fasta in options.fasta]
            inputFastaNames = [os.path.basename(fasta) for fasta in options.fasta]

            inputVCFFileIDs = []
            inputVCFNames = []
            inputTBIFileIDs = []
            if options.vcf:
                for vcf in options.vcf:
                    inputVCFFileIDs.append(toil.importFile(vcf))
                    inputVCFNames.append(os.path.basename(vcf))
                    inputTBIFileIDs.append(toil.importFile(vcf + '.tbi'))

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)

            # Unzip the fasta
            for i, fasta in enumerate(options.fasta):
                if fasta.endswith('.gz'):
                    inputFastaFileIDs[i] = init_job.addChildJobFn(run_unzip_fasta, context, inputFastaFileIDs[i], 
                                                                  os.path.basename(fasta)).rv()
                    inputFastaNames[i] = inputFastaNames[i][:-3]

            # Extract fasta sequence names and append them to regions
            if options.fasta_regions:
                scrape_fasta_job = init_job.addFollowOnJobFn(run_scan_fasta_sequence_names, context,
                                                             inputFastaFileIDs[0],
                                                             inputFastaNames[0],
                                                             options.regions)
                cur_job = scrape_fasta_job
                regions = scrape_fasta_job.rv()
            else:
                cur_job = init_job
                regions = options.regions

            # Automatically make and name a bunch of vcfs
            vcf_job = cur_job.addFollowOnJobFn(run_generate_input_vcfs, context, options.control_sample,
                                               inputVCFFileIDs, inputVCFNames,
                                               inputTBIFileIDs, 
                                               regions,
                                               options.out_name,
                                               filter_samples,
                                               options.haplo_sample,
                                               options.primary,
                                               options.min_af,
                                               not options.no_base)
                
            # Cosntruct graphs
            vcf_job.addFollowOnJobFn(run_construct_all, context, inputFastaFileIDs,
                                     inputFastaNames, vcf_job.rv(),
                                     options.max_node_size, options.alt_paths,
                                     options.flat_alts, regions,
                                     merge_graphs = options.merge_graphs,
                                     sort_ids = True, join_ids = True,
                                     gcsa_index = options.gcsa_index or options.all_index,
                                     xg_index = options.xg_index or options.all_index,
                                     gbwt_index = options.gbwt_index or options.all_index,
                                     snarls_index = options.snarls_index or options.all_index,
                                     haplo_sample = options.haplo_sample,
                                     gbwt_prune = options.gbwt_prune)
            
            # Run the workflow
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
