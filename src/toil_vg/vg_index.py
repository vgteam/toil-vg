#!/usr/bin/env python2.7
"""
vg_index.py: index a graph so it can be mapped to

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

logger = logging.getLogger(__name__)

def index_subparser(parser):
    """
    Create a subparser for indexing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")

    parser.add_argument("--graphs", nargs='+', type=make_url, required=True,
                        help="input graph(s). one per chromosome (separated by space)")

    parser.add_argument("--chroms", nargs='+',
                        help="Name(s) of reference path in graph(s) (separated by space).  If --graphs "
                        " has multiple elements, must be same length/order as --chroms (not needed for xg_index)")

    parser.add_argument("--node_mapping", type=make_url,
                        help="node mapping file required for gbwt pruning.  created by toil-vg construct"
                        " (or vg ids -j)")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add indexing options
    index_toggle_parse_args(parser)
    index_parse_args(parser)
    
    # Add common docker options
    add_container_tool_parse_args(parser)

def index_toggle_parse_args(parser):
    """ common args we do not want to share with toil-vg run """
    parser.add_argument("--gcsa_index", action="store_true",
                        help="Make a gcsa index for each output graph")
    parser.add_argument("--xg_index", action="store_true",
                        help="Make an xg index for each output graph")
    parser.add_argument("--gbwt_index", action="store_true",
                        help="Make a GBWT index alongside the xg index for each output graph")
    parser.add_argument("--snarls_index", action="store_true",
                        help="Make an snarls file for each output graph")
    parser.add_argument("--id_ranges_index", action="store_true",
                        help="Make chromosome id ranges tables (so toil-vg map can optionally split output by chromosome)")
    parser.add_argument("--all_index", action="store_true",
                        help="Equivalent to --gcsa_index --xg_index --gbwt_index --snarls_index --id_ranges_index")
    
def index_parse_args(parser):
    """ centralize indexing parameters here """
    
    parser.add_argument("--gcsa_index_cores", type=int,
        help="number of threads during the gcsa indexing step")
    parser.add_argument("--xg_index_cores", type=int,
        help="number of threads during the xg indexing step")
    parser.add_argument("--gbwt_index_cores", type=int,
        help="number of threads during the gbwt indexing step")    

    parser.add_argument("--kmers_cores", type=int,
        help="number of threads during the gcsa kmers step")

    parser.add_argument("--index_name", type=str, default='index',
                        help="name of index files. <name>.xg, <name>.gcsa etc.")

    parser.add_argument("--kmers_opts", type=str,
                        help="Options to pass to vg kmers.")
    parser.add_argument("--gcsa_opts", type=str,
                        help="Options to pass to gcsa indexing.")

    parser.add_argument("--vcf_phasing", nargs='+', type=make_url, default=[],
                        help="Import phasing information from VCF(s) into xg (or GBWT with --gbwt_index)")
    parser.add_argument("--gbwt_input", type=make_url,
                        help="Use given GBWT for GCSA2 pruning")
    parser.add_argument("--gbwt_prune", action='store_true',
                        help="Use gbwt for gcsa pruning")
                        
def validate_index_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    if any([options.gcsa_index, options.snarls_index,
            options.id_ranges_index, options.gbwt_index, options.all_index]):
        require(options.chroms, '--chroms must be specified')
        require(len(options.graphs) == 1 or len(options.chroms) == len(options.graphs),
                '--chroms and --graphs must have'
                ' same number of arguments if more than one graph specified if doing anything but xg indexing')
    require(any([options.xg_index, options.gcsa_index, options.snarls_index,
                 options.id_ranges_index, options.gbwt_index, options.all_index]),
            'at least one of --xg_index, --gcsa_index, --snarls_index, --id_ranged_index, --gbwt_index required, --all_index')
    require(not options.gbwt_prune or options.node_mapping,
                '--node_mapping required with --gbwt_prune')
    if options.vcf_phasing:
        require(all([vcf.endswith('.vcf.gz') for vcf in options.vcf_phasing]),
                'input phasing files must end with .vcf.gz')
    if options.gbwt_index:
        require(options.vcf_phasing, 'generating a GBWT requires a VCF with phasing information')
    if options.gbwt_prune:
        require(options.gbwt_index or options.gbwt_input, '--gbwt_index or --gbwt_input required for --gbwt_prune')
        require(options.gcsa_index or options.all_index, '--gbwt_prune requires gbwt indexing')
    require(not options.gbwt_index or not options.gbwt_input,
            'only one of --gbwt_index and --gbwt_input can be used at a time')
    if options.gbwt_input:
        require(options.gbwt_prune == 'gbwt', '--gbwt_prune required with --gbwt_input')
    
def run_gcsa_prune(job, context, graph_name, input_graph_id, gbwt_id, mapping_id):
    """
    Make a pruned graph using vg prune.  If unfold_mapping_id is provided, use -u, else -r
    """
    RealtimeLogger.info("Starting GBWT graph-pruning for gcsa kmers...")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Intermediate output
    restored_filename = os.path.join(work_dir, "restored_{}".format(graph_name))
    # Final output
    pruned_filename = os.path.join(work_dir, "unfolded_{}".format(graph_name))
    # Node Mapping output
    mapping_filename = os.path.join(work_dir, 'node_mapping')

    # Download input 
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(input_graph_id, graph_filename)
    gbwt_filename = graph_filename + '.gbwt'
    if gbwt_id:
        job.fileStore.readGlobalFile(gbwt_id, gbwt_filename)    
    if mapping_id:
        job.fileStore.readGlobalFile(mapping_id, mapping_filename, mutable=True)

    cmd = ['vg', 'prune', os.path.basename(graph_filename), '--threads', str(job.cores)]
    if context.config.prune_opts:
        cmd += context.config.prune_opts
    if gbwt_id:
        cmd += ['--append-mapping', '--mapping', os.path.basename(mapping_filename), '--unfold-paths']
        cmd += ['--gbwt-name', os.path.basename(gbwt_filename)]
    else:
        cmd += ['--restore-paths']
        
    with open(pruned_filename, 'w') as pruned_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile=pruned_file)
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GBWT pruning. Process took {} seconds.".format(run_time))

    pruned_graph_id = context.write_intermediate_file(job, pruned_filename)
    if gbwt_id:
        mapping_id = context.write_intermediate_file(job, mapping_filename)
    else:
        mapping_id = None

    return pruned_graph_id, mapping_id

def run_gcsa_kmers(job, context, graph_name, input_graph_id):
    """
    Make the kmers file, return its id
    """
    RealtimeLogger.info("Starting gcsa kmers...")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download input graph
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(input_graph_id, graph_filename)

    # Output
    output_kmers_filename = graph_filename + '.kmers'
   
    RealtimeLogger.info("Finding kmers in {} to {}".format(
        graph_filename, output_kmers_filename))

    # Make the GCSA2 kmers file
    with open(output_kmers_filename, "w") as kmers_file:
        command = ['vg', 'kmers',  os.path.basename(graph_filename), '--threads', str(job.cores)]
        command += context.config.kmers_opts
        context.runner.call(job, command, work_dir=work_dir, outfile=kmers_file)

    # Back to store
    output_kmers_id = context.write_intermediate_file(job, output_kmers_filename)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GCSA kmers. Process took {} seconds.".format(run_time))

    return output_kmers_id

def run_gcsa_prep(job, context, input_graph_ids,
                  graph_names, index_name, chroms,
                  chrom_gbwt_ids, node_mapping_id, skip_pruning=False):
    """
    Do all the preprocessing for gcsa indexing (pruning and kmers)
    Then launch the indexing as follow-on
    """    
    RealtimeLogger.info("Starting gcsa preprocessing...")
    start_time = timeit.default_timer()
    if chrom_gbwt_ids:
        assert len(chrom_gbwt_ids) <= len(input_graph_ids)

    kmers_ids = []

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)
    prune_root_job = Job()
    child_job.addChild(prune_root_job)

    # keep these in lists for now just in case
    prune_jobs = []
    # todo: figure out how best to update file with toil without making copies
    mapping_ids = [node_mapping_id] if node_mapping_id else []

    # prune then do kmers of each input graph.
    for graph_i, input_graph_id in enumerate(input_graph_ids):
        gbwt_id = chrom_gbwt_ids[graph_i] if chrom_gbwt_ids else None
        mapping_id = mapping_ids[-1] if mapping_ids else None        
        prev_job = prune_jobs[-1] if prune_jobs and gbwt_id else prune_root_job
        if not skip_pruning:
            prune_job = prev_job.addFollowOnJobFn(run_gcsa_prune, context, graph_names[graph_i],
                                                  input_graph_id, gbwt_id, mapping_id,
                                                  cores=context.config.prune_cores,
                                                  memory=context.config.prune_mem,
                                                  disk=context.config.prune_disk)
            prune_id = prune_job.rv(0)
            # toggle between parallele/sequential based on if we're unfolding or now
            if gbwt_id:
                prune_jobs.append(prune_job)
                mapping_ids.append(prune_job.rv(1))
        else:
            prune_id = input_graph_id

        # Compute the kmers as a follow-on to prune
        kmers_id = prune_job.addFollowOnJobFn(run_gcsa_kmers, context, graph_names[graph_i],
                                              prune_id, 
                                              cores=context.config.kmers_cores,
                                              memory=context.config.kmers_mem,
                                              disk=context.config.kmers_disk).rv()
        kmers_ids.append(kmers_id)

    return child_job.addFollowOnJobFn(run_gcsa_indexing, context, kmers_ids,
                                      graph_names, index_name, mapping_ids[-1] if mapping_ids else None,
                                      cores=context.config.gcsa_index_cores,
                                      memory=context.config.gcsa_index_mem,
                                      disk=context.config.gcsa_index_disk).rv()
    
def run_gcsa_indexing(job, context, kmers_ids, graph_names, index_name, mapping_id):
    """
    Make the gcsa2 index. Return its store id
    """
    
    RealtimeLogger.info("Starting gcsa indexing...")
    start_time = timeit.default_timer()     

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download all the kmers.  
    kmers_filenames = []
    
    for graph_i, kmers_id in enumerate(kmers_ids):
        kmers_filename = os.path.join(work_dir, os.path.basename(graph_names[graph_i]) + '.kmers')
        job.fileStore.readGlobalFile(kmers_id, kmers_filename)
        kmers_filenames.append(kmers_filename)

    # Download the mapping_id
    mapping_filename = None
    if mapping_id:
        mapping_filename = os.path.join(work_dir, 'node_mapping')
        job.fileStore.readGlobalFile(mapping_id, mapping_filename)

    # Where do we put the GCSA2 index?
    gcsa_filename = "{}.gcsa".format(index_name)

    command = ['vg', 'index', '--gcsa-name', os.path.basename(gcsa_filename)] + context.config.gcsa_opts
    command += ['--threads', str(job.cores)]
    
    for kmers_filename in kmers_filenames:
        command += ['--dbg-in', os.path.basename(kmers_filename)]
    if mapping_id:
        command += ['--mapping', os.path.basename(mapping_filename)]
    context.runner.call(job, command, work_dir=work_dir)

    # Checkpoint index to output store
    gcsa_file_id = context.write_output_file(job, os.path.join(work_dir, gcsa_filename))
    lcp_file_id = context.write_output_file(job, os.path.join(work_dir, gcsa_filename) + ".lcp")

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GCSA index. Process took {} seconds.".format(run_time))

    return gcsa_file_id, lcp_file_id

def run_concat_vcfs(job, context, vcf_ids, tbi_ids):
    """
    concatenate a list of vcfs.  we do this because vg index -v only takes one vcf, and
    we may be working with a set of chromosome vcfs. 
    """

    work_dir = job.fileStore.getLocalTempDir()

    vcf_names = ['chrom_{}.vcf.gz'.format(i) for i in range(len(vcf_ids))]
    out_name = 'genome.vcf.gz'

    for vcf_id, tbi_id, vcf_name in zip(vcf_ids, tbi_ids, vcf_names):
        job.fileStore.readGlobalFile(vcf_id, os.path.join(work_dir, vcf_name))
        job.fileStore.readGlobalFile(tbi_id, os.path.join(work_dir, vcf_name + '.tbi'))

    cmd = ['bcftools', 'concat'] + [vcf_name for vcf_name in vcf_names] + ['-O', 'z']
    with open(os.path.join(work_dir, out_name), 'w') as out_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile = out_file)

    cmd = ['tabix', '-f', '-p', 'vcf', out_name]
    context.runner.call(job, cmd, work_dir=work_dir)

    out_vcf_id = context.write_intermediate_file(job, os.path.join(work_dir, out_name))
    out_tbi_id = context.write_intermediate_file(job, os.path.join(work_dir, out_name + '.tbi'))

    return out_vcf_id, out_tbi_id

def run_xg_indexing(job, context, inputGraphFileIDs, graph_names, index_name,
                    vcf_phasing_file_id = None, tbi_phasing_file_id = None, make_gbwt = False):
    """
    Make the xg index and optional GBWT haplotype index.
    
    Saves the xg in the outstore as <index_name>.xg and the GBWT, if requested, as <index_name>.gbwt.
    
    Return a pair of file IDs, (xg_id, gbwt_id). The GBWT ID will be None if no GBWT is generated.
    """
    
    RealtimeLogger.info("Starting xg indexing...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Our local copy of the graphs
    graph_filenames = []
    for i, graph_id in enumerate(inputGraphFileIDs):
        graph_filename = os.path.join(work_dir, graph_names[i])
        job.fileStore.readGlobalFile(graph_id, graph_filename)
        graph_filenames.append(os.path.basename(graph_filename))

    # If we have a separate GBWT it will go here
    gbwt_filename = os.path.join(work_dir, "{}.gbwt".format(index_name))

    # Get the vcf file for loading phasing info
    if vcf_phasing_file_id:
        phasing_file = os.path.join(work_dir, 'phasing.vcf.gz')
        job.fileStore.readGlobalFile(vcf_phasing_file_id, phasing_file)
        job.fileStore.readGlobalFile(tbi_phasing_file_id, phasing_file + '.tbi')
        phasing_opts = ['-v', os.path.basename(phasing_file)]
        
        if make_gbwt and vcf_phasing_file_id:
            # Write the haplotype index to its own file
            phasing_opts += ['--gbwt-name', os.path.basename(gbwt_filename)]
    else:
        phasing_opts = []

    # Where do we put the XG index?
    xg_filename = "{}.xg".format(index_name)

    # Now run the indexer.
    RealtimeLogger.info("XG Indexing {}".format(str(graph_filenames)))

    command = ['vg', 'index', '--threads', str(job.cores), '--xg-name', os.path.basename(xg_filename)]
    command += phasing_opts + graph_filenames
    
    try:
        context.runner.call(job, command, work_dir=work_dir)
    except:
        # Dump everything we need to replicate the index run
        logging.error("XG indexing failed. Dumping files.")

        for graph_filename in graph_filenames:
            context.write_output_file(job, os.path.join(work_dir, graph_filename))
        if vcf_phasing_file_id:
            context.write_output_file(job, phasing_file)
            context.write_output_file(job, phasing_file + '.tbi')

        raise

    # Checkpoint index to output store
    xg_file_id = context.write_output_file(job, os.path.join(work_dir, xg_filename))
    
    gbwt_file_id = None
    if make_gbwt and vcf_phasing_file_id:
        # Also save the GBWT if it was generated
        gbwt_file_id = context.write_output_file(job, gbwt_filename)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished XG index. Process took {} seconds.".format(run_time))

    return (xg_file_id, gbwt_file_id)
    
def run_snarl_indexing(job, context, inputGraphFileIDs, graph_names, index_name):
    """
    Compute the snarls of the graph.
    
    Saves the snarls file in the outstore as <index_name>.snarls.
    
    Return the file ID of the snarls file.
    """
    
    RealtimeLogger.info("Starting snarl computation...")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Our local copy of the graphs
    graph_filenames = []
    for i, graph_id in enumerate(inputGraphFileIDs):
        graph_filename = os.path.join(work_dir, graph_names[i])
        job.fileStore.readGlobalFile(graph_id, graph_filename)
        graph_filenames.append(os.path.basename(graph_filename))

    # Where do we put the snarls?
    snarl_filename = os.path.join(work_dir, "{}.snarls".format(index_name))

    # Now run the indexer.
    RealtimeLogger.info("Computing Snarls for {}".format(str(graph_filenames)))

    pipeline = [['cat'] + graph_filenames, ['vg', 'snarls', '-']]
   
    with open(snarl_filename, "w") as snarl_file:
        # Concatenate all the graphs and compute the snarls.
        # Make sure to do it all in the container for vg (and not for 'cat')
        context.runner.call(job, pipeline, work_dir=work_dir, tool_name='vg', outfile=snarl_file)

    # Checkpoint index to output store
    snarl_file_id = context.write_output_file(job, snarl_filename)
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished Computing Snarls. Process took {} seconds.".format(run_time))

    return snarl_file_id


def run_id_ranges(job, context, inputGraphFileIDs, graph_names, index_name, chroms):
    """ Make a file of chrom_name <tab> first_id <tab> last_id covering the 
    id ranges of all chromosomes.  This is to speed up gam splitting down the road. 
    """
    
    RealtimeLogger.info("Starting id ranges...")
    start_time = timeit.default_timer()
    
    # Our id ranges (list of triples)
    id_ranges = []

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)    

    # Get the range for one graph per job. 
    for graph_id, graph_name, chrom in zip(inputGraphFileIDs, graph_names, chroms):
        id_range = child_job.addChildJobFn(run_id_range, context, graph_id, graph_name, chrom,
                                           cores=context.config.prune_cores,
                                           memory=context.config.prune_mem, disk=context.config.prune_disk).rv()
        
        id_ranges.append(id_range)

    # Merge them into a file and return its id
    return child_job.addFollowOnJobFn(run_merge_id_ranges, context, id_ranges, index_name,
                                      cores=context.config.misc_cores, memory=context.config.misc_mem,
                                      disk=context.config.misc_disk).rv()

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished id ranges. Process took {} seconds.".format(run_time))
    
def run_id_range(job, context, graph_id, graph_name, chrom):
    """
    Compute a node id range for a graph (which should be an entire contig/chromosome with
    contiguous id space -- see vg ids) using vg stats
    """
    work_dir = job.fileStore.getLocalTempDir()

    # download graph
    graph_filename = os.path.join(work_dir, graph_name)
    job.fileStore.readGlobalFile(graph_id, graph_filename)

    #run vg stats
    #expect result of form node-id-range <tab> first:last
    command = ['vg', 'stats', '--node-id-range', os.path.basename(graph_filename)]
    stats_out = context.runner.call(job, command, work_dir=work_dir, check_output = True).strip().split()
    assert stats_out[0] == 'node-id-range'
    first, last = stats_out[1].split(':')

    return chrom, first, last
    
def run_merge_id_ranges(job, context, id_ranges, index_name):
    """ create a BED-style file of id ranges
    """
    work_dir = job.fileStore.getLocalTempDir()

    # Where do we put the id ranges tsv?
    id_range_filename = os.path.join(work_dir, '{}_id_ranges.tsv'.format(index_name))

    with open(id_range_filename, 'w') as f:
        for id_range in id_ranges:
            f.write('{}\t{}\t{}\n'.format(*id_range))

    # Checkpoint index to output store
    return context.write_output_file(job, id_range_filename)

def run_merge_gbwts(job, context, chrom_gbwt_ids, index_name):
    """ merge up some gbwts
    """
    work_dir = job.fileStore.getLocalTempDir()

    gbwt_chrom_filenames = []

    for i, gbwt_id in enumerate(chrom_gbwt_ids):
        if gbwt_id:
            gbwt_filename = os.path.join(work_dir, '{}.gbwt'.format(i))
            job.fileStore.readGlobalFile(gbwt_id, gbwt_filename)
            gbwt_chrom_filenames.append(gbwt_filename)

    if len(gbwt_chrom_filenames) == 0:
        return None
    elif len(gbwt_chrom_filenames) == 1:
        return context.write_output_file(job, gbwt_chrom_filenames[0],
                                         out_store_path = index_name + '.gbwt')
    else:
        cmd = ['vg', 'gbwt', '--merge', '--fast', '--output', index_name + '.gbwt']
        cmd += [os.path.basename(f) for f in gbwt_chrom_filenames]
        context.runner.call(job, cmd, work_dir=work_dir)
        return context.write_output_file(job, os.path.join(work_dir, index_name + '.gbwt'))

def run_indexing(job, context, inputGraphFileIDs,
                 graph_names, index_name, chroms,
                 vcf_phasing_file_ids = [], tbi_phasing_file_ids = [], gbwt_id = None,
                 node_mapping_id = None,
                 skip_xg=False, skip_gcsa=False, skip_id_ranges=False,
                 skip_snarls=False, make_gbwt=False, gbwt_prune=False):
    """
    Run indexing logic by itself.
    
    Return a dict from index type ('xg','chrom_xg', 'gcsa', 'lcp', 'gbwt',
    'chrom_gbwt', 'id_ranges', or 'snarls') to index file ID(s) if created.
    
    For 'chrom_xg' and 'chrom_gbwt', the value is a list of one XG or GBWT per
    chromosome in chroms, to support `vg prune`. The others are all single file
    IDs
    
    """
    child_job = Job()
    job.addChild(child_job)
    xg_root_job = Job()
    child_job.addChild(xg_root_job)
    chrom_xg_root_job = Job()
    child_job.addChild(chrom_xg_root_job)

    # This will hold the index to return
    indexes = {}
    if gbwt_id:
        indexes['gbwt'] = gbwt_id

    make_gpbwt = vcf_phasing_file_ids and not make_gbwt
    
    if not skip_xg or not skip_gcsa:
        indexes['chrom_xg'] = []
        indexes['chrom_gbwt'] = []                                                            
        if make_gbwt:
            # In its current state, vg prune requires chromosomal xgs, so we must make
            # these xgs if we're doing any kind of gcsa indexing.  Also, if we're making
            # a gbwt, we do that at the same time (merging later if more than one graph)
            if not chroms or len(chroms) == 1:
                chroms = [index_name]
            indexes['chrom_xg'] = []
            indexes['chrom_gbwt'] = []                                        
            for i, chrom in enumerate(chroms):
                vcf_id = vcf_phasing_file_ids[i] if i < len(vcf_phasing_file_ids) else None
                tbi_id = tbi_phasing_file_ids[i] if i < len(tbi_phasing_file_ids) else None
                xg_chrom_index_job = chrom_xg_root_job.addChildJobFn(run_xg_indexing,
                                                                     context, [inputGraphFileIDs[i]],
                                                                     [graph_names[i]], chrom,
                                                                     vcf_id, tbi_id,
                                                                     make_gbwt = make_gbwt,
                                                                     cores=context.config.gbwt_index_cores,
                                                                     memory=context.config.gbwt_index_mem,
                                                                     disk=context.config.gbwt_index_disk)
                indexes['chrom_xg'].append(xg_chrom_index_job.rv(0))
                indexes['chrom_gbwt'].append(xg_chrom_index_job.rv(1))

            if len(chroms) > 1 and vcf_phasing_file_ids and make_gbwt:
                indexes['gbwt'] = chrom_xg_root_job.addFollowOnJobFn(run_merge_gbwts, context, indexes['chrom_gbwt'],
                                                                     index_name,
                                                                     cores=context.config.xg_index_cores,
                                                                     memory=context.config.xg_index_mem,
                                                                     disk=context.config.xg_index_disk).rv()

        # now do the whole genome xg (without any gbwt)
        if indexes.has_key('chrom_xg') and len(indexes['chrom_xg']) == 1 and not make_gpbwt:
            # our first chromosome is effectively the whole genome (note that above we
            # detected this and put in index_name so it's saved right (don't care about chrom names))
            indexes['xg'] = indexes['chrom_xg'][0]
        else:            
            if make_gpbwt and len(vcf_phasing_file_ids) > 1:
                concat_job = xg_root_job.addChildJobFn(run_concat_vcfs, context,
                                                       vcf_phasing_file_ids, tbi_phasing_file_ids,
                                                       cores=1,
                                                       memory=context.config.xg_index_mem,
                                                       disk=context.config.xg_index_disk)
                vcf_phasing_file_id = concat_job.rv(0)
                tbi_phasing_file_id = concat_job.rv(1)
            else:
                concat_job = Job()
                xg_root_job.addChild(concat_job)
                vcf_phasing_file_id = None
                tbi_phasing_file_id = None

            if not skip_xg:
                xg_index_job = concat_job.addChildJobFn(run_xg_indexing,
                                                        context, inputGraphFileIDs,
                                                        graph_names, index_name,
                                                        vcf_phasing_file_id, tbi_phasing_file_id,
                                                        make_gbwt=False,
                                                        cores=context.config.xg_index_cores,
                                                        memory=context.config.xg_index_mem,
                                                        disk=context.config.xg_index_disk)
                indexes['xg'] = xg_index_job.rv(0)


    gcsa_root_job = Job()
    # gcsa follow from chrom_xg jobs only if gbwt needed for pruning
    gcsa_predecessor_job = chrom_xg_root_job if gbwt_prune else child_job
    gcsa_predecessor_job.addFollowOn(gcsa_root_job)
    
    if not skip_gcsa:
        # We know we made the per-chromosome indexes already, so we can use them here to make the GCSA                                               
        # todo: we're only taking in a genome gbwt as input, because that's all we write
        if not indexes.has_key('chrom_gbwt') and indexes.has_key('gbwt'):
            indexes['chrom_gbwt'] = indexes['gbwt'] * len(cinputGraphFileIDs)
        gcsa_job = gcsa_root_job.addChildJobFn(run_gcsa_prep, context, inputGraphFileIDs,
                                               graph_names, index_name, chroms,
                                               indexes['chrom_gbwt'] if gbwt_prune else [],
                                               node_mapping_id,
                                               cores=context.config.misc_cores,
                                               memory=context.config.misc_mem,
                                               disk=context.config.misc_disk)
        indexes['gcsa'] = gcsa_job.rv(0)
        indexes['lcp'] = gcsa_job.rv(1)
    
    if len(inputGraphFileIDs) > 1 and not skip_id_ranges:
        indexes['id_ranges'] = job.addChildJobFn(run_id_ranges, context, inputGraphFileIDs,
                                                 graph_names, index_name, chroms,
                                                 cores=context.config.misc_cores,
                                                 memory=context.config.misc_mem,
                                                 disk=context.config.misc_disk).rv()

    return indexes

def index_main(context, options):
    """
    Wrapper for vg indexing. 
    """

    # check some options
    validate_index_options(options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            inputGraphFileIDs = []
            for graph in options.graphs:
                inputGraphFileIDs.append(toil.importFile(graph))
            inputPhasingVCFFileIDs = []
            inputPhasingTBIFileIDs = []
            for vcf in options.vcf_phasing:
                inputPhasingVCFFileIDs.append(toil.importFile(vcf))
                inputPhasingTBIFileIDs.append(toil.importFile(vcf + '.tbi'))
            # like construct, if there's one vcf and many graphs, we apply the vcf everywhere
            if len(inputPhasingTBIFileIDs) == 1 and len(options.graphs) > 1:
                inputPhasingVCFFileIDs = [inputPhasingVCFFileIDs[0]] * len(options.graphs)
                inputPhasingTBIFileIDs = [inputPhasingTBIFileIDs[0]] * len(options.graphs)
            inputGBWTID = None
            if options.gbwt_input:
                inputGBWTID = toil.importFile(options.gbwt_input)
            inputNodeMappingID = None
            if options.node_mapping:
                inputNodeMappingID = toil.importFile(options.node_mapping)

            # Handy to have meaningful filenames throughout, so we remember
            # the input graph names
            graph_names = [os.path.basename(i) for i in options.graphs]

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_indexing, context, inputGraphFileIDs,
                                     graph_names, options.index_name, options.chroms,
                                     inputPhasingVCFFileIDs, inputPhasingTBIFileIDs,
                                     inputGBWTID, inputNodeMappingID,
                                     skip_xg = not options.xg_index and not options.all_index,
                                     skip_gcsa = not options.gcsa_index and not options.all_index,
                                     skip_id_ranges = not options.id_ranges_index and not options.all_index,
                                     skip_snarls = not options.snarls_index and not options.all_index,
                                     make_gbwt=options.gbwt_index, gbwt_prune=options.gbwt_prune,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)
            init_job.addFollowOn(root_job)            
            
            # Run the job and store the returned list of output files to download
            index_key_and_id = toil.start(init_job)
        else:
            index_key_and_id = toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
