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

    parser.add_argument("--skip_xg", action="store_true",
                        help="Do not generate xg index")
    parser.add_argument("--skip_gcsa", action="store_true",
                        help="Do not generate gcsa index")
    parser.add_argument("--skip_id_ranges", action="store_true",
                        help="Do not generate id_ranges.tsv")
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add indexing options
    index_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)


def index_parse_args(parser):
    """ centralize indexing parameters here """

    parser.add_argument("--graphs", nargs='+', type=make_url,
                        help="input graph(s). one per chromosome (separated by space)")

    parser.add_argument("--chroms", nargs='+',
                        help="Name(s) of reference path in graph(s) (separated by space).  If --graphs "
                        " has multiple elements, must be same length/order as --chroms")

    parser.add_argument("--gcsa_index_cores", type=int,
        help="number of threads during the gcsa indexing step")
    parser.add_argument("--xg_index_cores", type=int,
        help="number of threads during the xg indexing step")

    parser.add_argument("--kmers_cores", type=int,
        help="number of threads during the gcsa kmers step")

    parser.add_argument("--index_name", type=str, default='index',
                        help="name of index files. <name>.xg, <name>.gcsa etc.")

    parser.add_argument("--kmers_opts", type=str,
                        help="Options to pass to vg kmers.")
    parser.add_argument("--gcsa_opts", type=str,
                        help="Options to pass to gcsa indexing.")

    parser.add_argument("--vcf_phasing", nargs='+', type=make_url, default=[],
                        help="Import phasing information from VCF(s) into xg or GBWT")
    parser.add_argument("--make_gbwt", action='store_true',
                        help="Save phasing information to a GBWT (instead of GBWT inside XG)")
    # todo: do we want an option to pass in a GBWT?
    parser.add_argument("--haplo_pruning", action='store_true',
                        help="Use GBWT for haplotype pruning for GCSA construction")
                        
def validate_index_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """                           
    require(options.chroms and options.graphs, '--chroms and --graphs must be specified')
    require(len(options.graphs) == 1 or len(options.chroms) == len(options.graphs),
            '--chroms and --graphs must have'
            ' same number of arguments if more than one graph specified')
    if options.vcf_phasing:
        require(all([vcf.endswith('.vcf.gz') for vcf in options.vcf_phasing]),
                'input phasing files must end with .vcf.gz')
    if options.make_gbwt:
        require(options.vcf_phasing, 'generating a GBWT requires a VCF with phasing information')
    if options.haplo_pruning:
        require(options.make_gbwt, '--make_gbwt required for --haplo_pruning')
    
def run_gcsa_prune(job, context, graph_name, input_graph_id, primary_paths=[]):
    """
    Make the pruned graph, return its id
    Retains the specified primary paths.
    """
    RealtimeLogger.info("Starting graph-pruning for gcsa kmers...")
    start_time = timeit.default_timer()

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Place where we put pruned vg
    to_index_filename = None
    
    if len(context.config.prune_opts) > 0:

        # Make sure we are dealing with list of lists
        if not isinstance(context.config.prune_opts[0], list):
            context.config.prune_opts = [context.config.prune_opts]
        assert all([isinstance(x, list) for x in  context.config.prune_opts])

        # Download input graph
        graph_filename = os.path.join(work_dir, graph_name)
        job.fileStore.readGlobalFile(input_graph_id, graph_filename)

        to_index_filename = os.path.join(work_dir, "pruned_{}".format(graph_name))
        with open(to_index_filename, "w") as to_index_file:
            command = [['vg', 'mod', os.path.basename(graph_filename), '-t',
                        str(max(1, job.cores / len(context.config.prune_opts)))] +
                       context.config.prune_opts[0]]
            
            # tack on subsequent vg mod commands if specified
            for po in context.config.prune_opts[1:]:
                command.append(['vg', 'mod', '-', '-t', 
                                str(max(1, job.cores / len(context.config.prune_opts)))] + po)
                               
            context.runner.call(job, command, work_dir=work_dir, outfile=to_index_file)
            
            # Then append in the primary path.
            command = ['vg', 'mod', '-N']
            for primary_path in primary_paths:
                # Send along all the primary paths
                command.append('-r')
                command.append(primary_path)
            command += ['-t', str(job.cores), os.path.basename(graph_filename)]
            context.runner.call(job, command, work_dir=work_dir, outfile=to_index_file)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished pruning. Process took {} seconds.".format(run_time))

    if to_index_filename is None:
        # no pruning done: just pass along the input graph as is
        pruned_graph_id = input_graph_id
    else:
        pruned_graph_id = context.write_intermediate_file(job, to_index_filename)

    return pruned_graph_id

def run_gcsa_prune_gbwt(job, context, graph_name, input_graph_id, xg_id, gbwt_id, mapping_id):
    """
    Make the pruned graph using gbwt and vg prune and return (pruned_graph_id, mapping_id)
    (Can be used even if gbwt_id is None)
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
    xg_filename = graph_filename + '.xg'
    job.fileStore.readGlobalFile(xg_id, xg_filename)
    gbwt_filename = graph_filename + '.gbwt'
    if gbwt_id:
        job.fileStore.readGlobalFile(gbwt_id, gbwt_filename)    
    if mapping_id:
        job.fileStore.readGlobalFile(mapping_id, mapping_filename, mutable=True)

    # (from https://github.com/jltsiren/gcsa2/wiki/New-Construction-Benchmarks)
    # first pruning step (default edge max restoring paths) TODO: parameterize
    cmd = ['vg', 'prune', '-r', '-x', os.path.basename(xg_filename),
           '-p', os.path.basename(graph_filename), '-t', str(job.cores)]
    with open(restored_filename, 'w') as restored_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile=restored_file)

    # second pruning step (remove unobserved recombinations using gbwt)
    cmd = ['vg', 'prune', '-u', '-x', os.path.basename(xg_filename),
           '-a',
           '-m', os.path.basename(mapping_filename),
           '-p', os.path.basename(graph_filename),
           '-t', str(job.cores)]
    if gbwt_id:
        cmd += ['-g', os.path.basename(gbwt_filename)]
    with open(pruned_filename, 'w') as pruned_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile=pruned_file)
    
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished GBWT pruning. Process took {} seconds.".format(run_time))

    pruned_graph_id = context.write_intermediate_file(job, pruned_filename)
    mapping_id = context.write_intermediate_file(job, mapping_filename)

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
        command = ['vg', 'kmers',  os.path.basename(graph_filename), '-t', str(job.cores)]
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
                  chrom_xg_ids=None, chrom_gbwt_ids=None,
                  primary_path_override=None):
    """
    Do all the preprocessing for gcsa indexing (pruning and kmers)
    Then launch the indexing as follow-on
    """    
    RealtimeLogger.info("Starting gcsa preprocessing...")
    start_time = timeit.default_timer()
    if chrom_gbwt_ids:
        assert len(chrom_gbwt_ids) == len(chrom_xg_ids)
        assert len(chrom_gbwt_ids) <= len(input_graph_ids)
    if chrom_xg_ids:
        assert len(chrom_xg_ids) == len(input_graph_ids)

    kmers_ids = []

    # to encapsulate everything under this job
    child_job = Job()
    job.addChild(child_job)
    prune_root_job = Job()
    child_job.addChild(prune_root_job)

    # keep these in lists for now just in case
    prune_jobs = []
    mapping_ids = []

    # prune then do kmers of each input graph.  note that new interface pruning needs to be
    # done sequentially
    for graph_i, input_graph_id in enumerate(input_graph_ids):
        # toggle new interface with chrom_xg_ids
        if chrom_xg_ids:
            # Do the new GBWT-based pruning
            gbwt_id = chrom_gbwt_ids[graph_i] if chrom_gbwt_ids else None
            xg_id = chrom_xg_ids[graph_i]
            prev_job = prune_jobs[-1] if prune_jobs else prune_root_job
            mapping_id = mapping_ids[-1] if mapping_ids else None
            prune_job = prev_job.addFollowOnJobFn(run_gcsa_prune_gbwt, context, graph_names[graph_i],
                                                  input_graph_id, xg_id, gbwt_id, mapping_id,
                                                  cores=context.config.prune_cores,
                                                  memory=context.config.prune_mem,
                                                  disk=context.config.prune_disk)
            prune_id = prune_job.rv(0)
            prune_jobs.append(prune_job)
            mapping_ids.append(prune_job.rv(1))

        else:
            # Determine the primary path list to use
            if primary_path_override is not None:
                primary_paths = primary_path_override
            elif len(input_graph_ids) == len(chroms):
                primary_paths = [chroms[graph_i]]
            else:
                assert len(input_graph_ids) == 1
                primary_paths = chroms

            # Do the old mod-based pruning where we keep a path
            # TODO: Can we try to use vg prune here too?
            prune_job = prune_root_job.addChildJobFn(run_gcsa_prune, context, graph_names[graph_i],
                                                     input_graph_id, primary_paths=primary_paths,
                                                     cores=context.config.prune_cores,
                                                     memory=context.config.prune_mem,
                                                     disk=context.config.prune_disk)
            prune_id = prune_job.rv()

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

    command = ['vg', 'index', '-g', os.path.basename(gcsa_filename)] + context.config.gcsa_opts
    command += ['-t', str(job.cores)]
    for kmers_filename in kmers_filenames:
        command += ['-i', os.path.basename(kmers_filename)]
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
        
        if make_gbwt:
            # Write the haplotype index to its own file
            phasing_opts += ['--gbwt-name', os.path.basename(gbwt_filename)]
    else:
        phasing_opts = []

    # Where do we put the XG index?
    xg_filename = "{}.xg".format(index_name)

    # Now run the indexer.
    RealtimeLogger.info("XG Indexing {}".format(str(graph_filenames)))

    command = ['vg', 'index', '-t', str(job.cores), '-x', os.path.basename(xg_filename)]
    command += phasing_opts + graph_filenames
    
    context.runner.call(job, command, work_dir=work_dir)

    # Checkpoint index to output store
    xg_file_id = context.write_output_file(job, os.path.join(work_dir, xg_filename))
    
    gbwt_file_id = None
    if make_gbwt:
        # Also save the GBWT if it was generated
        gbwt_file_id = context.write_output_file(job, gbwt_filename)

    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Finished XG index. Process took {} seconds.".format(run_time))

    return (xg_file_id, gbwt_file_id)


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
    command = ['vg', 'stats', '-r', os.path.basename(graph_filename)]
    stats_out = context.runner.call(job, command, work_dir=work_dir, check_output = True).strip().split()
    assert stats_out[0] == 'node-id-range'
    first, last = stats_out[1].split(':')

    return chrom, first, last
    
def run_merge_id_ranges(job, context, id_ranges, index_name):
    """ create a BED-style file of id ranges
    """
    work_dir = job.fileStore.getLocalTempDir()

    # Where do we put the XG index?
    id_range_filename = os.path.join(work_dir, '{}_id_ranges.tsv'.format(index_name))

    with open(id_range_filename, 'w') as f:
        for id_range in id_ranges:
            f.write('{}\t{}\t{}\n'.format(*id_range))

    # Checkpoint index to output store
    return context.write_output_file(job, id_range_filename)

def run_indexing(job, context, inputGraphFileIDs,
                 graph_names, index_name, chroms,
                 vcf_phasing_file_ids = [], tbi_phasing_file_ids = [],
                 skip_xg=False, skip_gcsa=False, skip_id_ranges=False, make_gbwt=False,
                 haplo_pruning=False):
    """
    Run indexing logic by itself.
    
    Return an XG file ID, a list of chrom XG IDs, a GBWT ID, a list of chrom GBWT IDs,
    a GCSA ID, a LCP ID, a GPBWT
    and an ID for the ID ranges index file. (any of the above can be empty depending on the flags)
    
    """
    xg_root_job = Job()
    job.addChild(xg_root_job)

    make_gpbwt = vcf_phasing_file_ids and not make_gbwt
    
    # return a tuple of this:
    xg_index_id = None
    chrom_xg_ids = []
    gbwt_id = None
    chrom_gbwt_ids = []
    gcsa_id = None
    lcp_id = None
    gpbwt_id = None
    id_ranges_id = None
    
    if not skip_xg:
        if make_gbwt:
            # need to make chrosomome xgs and gbwts.  there is annoying redundance in
            # xg indexing as we'll make a whole-genome xg, but there's no real getting
            # around it atm (can add options to do one or there other if necessary).
            if not chroms or len(chroms) == 1:
                chroms = [index_name]
            for i, chrom in enumerate(chroms):
                vcf_id = vcf_phasing_file_ids[i] if i < len(vcf_phasing_file_ids) else None
                tbi_id = tbi_phasing_file_ids[i] if i < len(tbi_phasing_file_ids) else None
                xg_chrom_index_job = xg_root_job.addChildJobFn(run_xg_indexing,
                                                               context, [inputGraphFileIDs[i]],
                                                               [graph_names[i]], chrom,
                                                               vcf_id, tbi_id,
                                                               make_gbwt = vcf_id is not None,
                                                               cores=context.config.xg_index_cores,
                                                               memory=context.config.xg_index_mem,
                                                               disk=context.config.xg_index_disk)
                chrom_xg_ids.append(xg_chrom_index_job.rv(0))
                chrom_gbwt_ids.append(xg_chrom_index_job.rv(1))

            if len(vcf_phasing_file_ids) == 1:
                gbwt_id = chrom_gbwt_ids[0]
            else:
                # todo: optional merge?
                gbwt_id = xg_root_job.addFollowOnJobFn(run_merge_gbwts, context, chrom_gbwt_ids,
                                                       cores=context.config.xg_index_cores,
                                                       memory=context.config.xg_index_mem,
                                                       disk=context.config.xg_index_disk).rv()

        # now do the whole genome xg (without any gbwt)
        if len(chrom_xg_ids) == 1 and not make_gpbwt:
            # our first chromosome is effectively the whole genome (note that above we
            # detected this and put in index_name so it's saved right (don't care about chrom names))
            xg_index_id = chrom_xg_ids[0]
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
                vcf_phasing_file_id = None if not vcf_phasing_file_ids else vcf_phasing_file_ids[0]
                tbi_phasing_file_id = None if not tbi_phasing_file_ids else tbi_phasing_file_ids[0]
                
            xg_index_job = concat_job.addChildJobFn(run_xg_indexing,
                                                    context, inputGraphFileIDs,
                                                    graph_names, index_name,
                                                    vcf_phasing_file_id, tbi_phasing_file_id,
                                                    make_gbwt=False,
                                                    cores=context.config.xg_index_cores,
                                                    memory=context.config.xg_index_mem,
                                                    disk=context.config.xg_index_disk)
            xg_index_id = xg_index_job.rv(0)

    # the haplo_pruning flag determines if gcsa is a follow-on from xg
    gcsa_root_job = Job()
    if haplo_pruning:
        xg_root_job.addFollowOn(gcsa_root_job)
    else:
        job.addChild(gcsa_root_job)
    
    if not skip_gcsa:
        gcsa_job = gcsa_root_job.addChildJobFn(run_gcsa_prep, context, inputGraphFileIDs,
                                               graph_names, index_name, chroms, chrom_xg_ids,
                                               chrom_gbwt_ids,
                                               cores=context.config.misc_cores,
                                               memory=context.config.misc_mem,
                                               disk=context.config.misc_disk)
        gcsa_id = gcsa_job.rv(0)
        lcp_id = gcsa_job.rv(1)
    
    if len(inputGraphFileIDs) > 1 and not skip_id_ranges:
        id_ranges_id = job.addChildJobFn(run_id_ranges, context, inputGraphFileIDs,
                                         graph_names, index_name, chroms,
                                         cores=context.config.misc_cores,
                                         memory=context.config.misc_mem,
                                         disk=context.config.misc_disk).rv()

    return xg_index_id, chrom_xg_ids, gbwt_id, chrom_gbwt_ids, gcsa_id, lcp_id, id_ranges_id


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

            # Handy to have meaningful filenames throughout, so we remember
            # the input graph names
            graph_names = [os.path.basename(i) for i in options.graphs]

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_indexing, context, inputGraphFileIDs,
                                     graph_names, options.index_name, options.chroms,
                                     inputPhasingVCFFileIDs, inputPhasingTBIFileIDs,
                                     options.skip_xg, options.skip_gcsa, options.skip_id_ranges,
                                     options.make_gbwt,
                                     options.haplo_pruning,
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
    
