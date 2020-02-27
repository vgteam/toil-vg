#!/usr/bin/env python
"""
vg_msga.py: use vg's banded aligner to align contigs into the graph.  By associating BED regions with each
contig, we can parallelize by chromosome. 

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

def msga_subparser(parser):
    """
    Create a subparser for msgaing.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--graphs", nargs='+', default=[], type=make_url,
                        help="input graph(s). one per chromosome (separated by space)")
    parser.add_argument("--chroms", nargs='+',
                        help="name(s) of reference path in graph(s) (separated by space).  If --graphs "
                        " has multiple elements, must be same length/order as --chroms (not needed for xg_index)")
    parser.add_argument("--target_regions", type=make_url,
                        help="BED file mapping regions (cols 1-3) to sequence names (col 4) from the FASTA")
    parser.add_argument("--fasta", type=make_url, required = True,
                        help="FASTA file containing sequences to align")

    # Add msga options
    msga_parse_args(parser)
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)

def msga_parse_args(parser):
    """
    Indexing parameters which can be part of construction pipeline
    """
    parser.add_argument("--alignment_cores", type=int,
                        help="number of threads during the alignment step")
    parser.add_argument("--msga_context", type=int,
                        help="number of context steps when expanding target region")

def validate_msga_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(len(options.graphs) == 0 or options.chroms, '--chroms must be specified for --graphs')
    require(len(options.chroms) == len(options.graphs),
            '--chroms and --graphs must have same number of arguments')
    require(len(options.chroms) <= 1 or options.target_regions,
             '--target_regions required for multiple chromosomes')

def run_msga(job, context, graph_name, graph_id, fasta_id, target_regions_id, chrom, normalize = False, max_node_size = 32, validate = False):
    """
    Run vg msga to align some fasta sequences to a graph
    """

    if chrom and isinstance(chrom, list):
        # we parallelize by recursing on each chromosome
        if not graph_id:
            graph_id = [None] * len(chroms)
            graph_path = [None] * len(chroms)
        return [job.addChildJobFn(run_msga, context, g_name, g_id, fasta_id, target_regions_id, g_chrom,
                                  cores=context.config.alignment_cores,
                                  memory=context.config.alignment_mem,
                                  disk=context.config.alignment_disk).rv() \
                for g_name, g_id, g_chrom in zip(graph_name, graph_id, chrom)]
    else:
        # run on one chromosome
        work_dir = job.fileStore.getLocalTempDir()

        # process our regions file
        target_regions = []
        if target_regions_id:
            regions_path = os.path.join(work_dir, 'regions.bed')
            job.fileStore.readGlobalFile(target_regions_id, regions_path)
            with open(regions_path) as regions_file:
                for line in regions_file:
                    toks = line.strip().split('\t')
                    if len(toks) >= 4 and not toks[0].startswith('#'):
                        bed_chrom, bed_name = toks[0], toks[3]
                        if not chrom or bed_chrom == chrom:
                            target_regions.append(bed_name)
            # Don't bother continuing if there's nothing to map                            
            if not target_regions:
                RealtimeLogger.info("No sequences found for chromosome {}".format(chrom))
                return graph_id

        # download input
        if graph_id:
            if isinstance(graph_id, list):
                assert len(graph_id) == 1 and len(graph_name) == 1
                graph_id = graph_id[0]
                graph_name = graph_name[0]
            graph_path = os.path.join(work_dir, graph_name)
            job.fileStore.readGlobalFile(graph_id, graph_path)
        else:
            graph_path = 'graph.vg'
        fasta_path = os.path.join(work_dir, 'contigs.fa')
        job.fileStore.readGlobalFile(fasta_id, fasta_path)

        # subset the fasta to only contain sequences that belong in chrom (as determined from the BED)
        if target_regions_id and chrom:
            fasta_subset_path = os.path.join(work_dir, '{}.fa'.format(chrom))
            subset_count = 0
            with open(fasta_subset_path, 'wb') as subset_file:
                for target_region in target_regions:
                    toks = line.strip().split('\t')
                    context.runner.call(job, ['samtools', 'faidx', os.path.basename(fasta_path),
                                              target_region], work_dir = work_dir, outfile = subset_file)
            fasta_path = fasta_subset_path

        # run vg msga to align each fasta sequence to our graph (iteratively, doing the longest first)
        msga_cmd = ['vg', 'msga', '--from', os.path.basename(fasta_path), '--threads', str(job.cores)]
        if graph_id:
            msga_cmd += ['--graph', os.path.basename(graph_path)]
        if target_regions_id:
            msga_cmd += ['--position-bed', os.path.basename(regions_path), '--context' ,
                         str(context.config.msga_context)]
        msga_cmd += context.config.msga_opts
        
        if normalize:
            # msga's built-in normalization is run incrementally which my improve alignment
            msga_cmd += ['--normalize']
            # but we still normalize at the end, in part just to make sure that the node size is
            # respected (only controlable in msga via kmer size parameter).  can't be too normal
            msga_cmd = [msga_cmd,
                        ['vg', 'mod', '--until-normal', str(context.config.normalize_iterations), '-']]
            msga_cmd.append(['vg', 'mod', '--chop', str(max_node_size), '-'])
            msga_cmd.append(['vg', 'ids', '--sort', '-'])
        
        out_path = graph_path[:-3] + '-msga.vg'
        try:
            with open(out_path, 'w') as out_file:
                context.runner.call(job, msga_cmd, work_dir = work_dir, outfile = out_file)
        except:
            # Dump everything we need to replicate the construction
            logging.error("msga failed. Dumping files.")
            for dump_path in [fasta_path, graph_path, regions_path]:
                if dump_path and os.path.isfile(dump_path):
                    context.write_output_file(job, dump_path)
            raise

        # Check the graph for errors
        if validate:
            context.runner.call(job, ['vg', 'validate', os.path.basename(out_path)], work_dir = work_dir)

        return context.write_output_file(job, out_path)
    
    
def msga_main(context, options):
    """
    Wrapper for vg msga. 
    """

    # check some options
    validate_msga_options(options)

    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Upload local files to the remote IO Store
            graph_ids = [importer.load(graph) for graph in options.graphs]
            target_regions_id = importer.load(options.target_regions) if options.target_regions else None
            fasta_id = importer.load(options.fasta) if options.fasta else None

            importer.wait()            
            
            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
            cur_job = init_job

            # Unzip the fasta
            if options.fasta.endswith('.gz'):
                fasta_id = cur_job.addChildJobFn(run_unzip_fasta, context, importer.resolve(fasta_id), 
                                                 os.path.basename(options.fasta),
                                                 disk=context.config.construct_disk).rv()
                #todo: name mangling and iupac fixing?

            # Make a root job
            cur_job = cur_job.addFollowOnJobFn(run_msga, context,
                                               [os.path.basename(graph) for graph in options.graphs],
                                               importer.resolve(graph_ids),
                                               importer.resolve(fasta_id),
                                               importer.resolve(target_regions_id),
                                               options.chroms,
                                               cores=context.config.alignment_cores,
                                               memory=context.config.alignment_mem,
                                               disk=context.config.alignment_disk)
            # Run the job and store the returned list of output files to download
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    logger.info("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
