#!/usr/bin/env python2.7
"""
vg_pedigree.py: pedigree map and calling pipeline to produce parental-enhanced mapping and calling output.

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
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
from toil_vg.vg_call import *
from toil_vg.vg_index import *
from toil_vg.vg_map import *
from toil_vg.vg_surject import *
from toil_vg.vg_config import *
from toil_vg.vg_construct import *
from toil_vg.context import Context, run_write_info_to_outstore

logger = logging.getLogger(__name__)

def pedigree_subparser(parser):
    """
    Create a subparser for pedigree workflow.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("proband_name", type=str,
                        help="sample name of proband or sample of interest (ex HG002)")
    parser.add_argument("maternal_name", type=str,
                        help="sample name of mother to the proband (ex HG004)")
    parser.add_argument("paternal_name", type=str,
                        help="sample name of father to the proband (ex HG003)")
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--sibling_names", nargs='+', type=str, default=None,
                        help="sample names of siblings to the proband. Optional.")
    parser.add_argument("--kmer_size", type=int,
                        help="size of kmers to use in gcsa-kmer mapping mode")
    parser.add_argument("--id_ranges", type=make_url, default=None,
                        help="Path to file with node id ranges for each chromosome in BED format.")
    parser.add_argument("--fastq_proband", nargs='+', type=make_url,
                        help="Proband input fastq(s) (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fastq_maternal", nargs='+', type=make_url,
                        help="Maternal input fastq(s) (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fastq_paternal", nargs='+', type=make_url,
                        help="Paternal input fastq(s) (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--fastq_siblings", nargs='+', type=make_url,
                        help="Sibling input fastq(s) (possibly compressed), two are allowed, one for each mate per sibling.
                            Sibling read-pairs must be input adjacent to eachother. Must follow same order as input to
                            --sibling_names argument.")
    parser.add_argument("--gam_input_reads_proband", type=make_url, default=None,
                        help="Input reads of proband in GAM format")
    parser.add_argument("--gam_input_reads_maternal", type=make_url, default=None,
                        help="Input reads of mother in GAM format")
    parser.add_argument("--gam_input_reads_paternal", type=make_url, default=None,
                        help="Input reads of father in GAM format")
    parser.add_argument("--gam_input_reads_siblings", nargs='+', type=make_url, default=None,
                        help="Input reads of sibling(s) in GAM format. Must follow same order as input to
                            --sibling_names argument.")
    parser.add_argument("--bam_input_reads_proband", type=make_url, default=None,
                        help="Input reads of proband in BAM format")
    parser.add_argument("--bam_input_reads_maternal", type=make_url, default=None,
                        help="Input reads of mother in BAM format")
    parser.add_argument("--bam_input_reads_paternal", type=make_url, default=None,
                        help="Input reads of father in BAM format")
    parser.add_argument("--bam_input_reads_siblings", nargs='+', type=make_url, default=None,
                        help="Input reads of sibling(s) in BAM format. Must follow same order as input to
                            --sibling_names argument.")
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add mapping index options
    map_parse_index_args(parser)

    # Add pedigree options shared only with map
    pedigree_parse_index_args(parser)
    
    # Add common docker options
    add_container_tool_parse_args(parser)

def pedigree_parse_args(parser, stand_alone = False):
    """
    Define pedigree arguments shared with map
    """

    parser.add_argument("--fq_split_cores", type=int,
                        help="number of threads used to split input FASTQs")
    parser.add_argument("--single_reads_chunk", action="store_true", default=False,
                        help="do not split reads into chunks")
    parser.add_argument("--reads_per_chunk", type=int,
                        help="number of reads for each mapping job")
    parser.add_argument("--alignment_cores", type=int,
                        help="number of threads during the alignment step")
    parser.add_argument("--interleaved", action="store_true", default=False,
                        help="treat fastq as interleaved read pairs. overrides *_opts")
    parser.add_argument("--map_opts", type=str,
                        help="arguments for vg map (wrapped in \"\")")
    parser.add_argument("--mpmap_opts", type=str,
                        help="arguments for vg mpmap (wrapped in \"\")")
    parser.add_argument("--gaffe_opts", type=str,
                        help="arguments for vg gaffe (wrapped in \"\")")
    parser.add_argument("--bam_output", action="store_true",
                        help="write BAM output directly")
    parser.add_argument("--surject", action="store_true",
                        help="surject output, producing BAM in addition to GAM alignments")
    parser.add_argument("--validate", action="store_true",
                        help="run vg validate on ouput GAMs")

def validate_pedigree_options(context, options):
    """
    Throw an error if an invalid combination of options has been selected.
    """
    require(options.xg_index is not None, 'All mappers require --xg_index')
    
    if options.mapper == 'map' or options.mapper == 'mpmap':
        require(options.gcsa_index, '--gcsa_index is required for map and mpmap')
    
    if options.mapper == 'gaffe':
        require(options.minimizer_index, '--minimizer_index is required for gaffe')
        require(options.distance_index, '--distance_index is required for gaffe')
        require(options.gbwt_index, '--gbwt_index is required for gaffe')
        require(not options.bam_input_reads, '--bam_input_reads is not supported with gaffe')
        require(not options.interleaved, '--interleaved is not supported with gaffe')
        require(options.fastq is None or len(options.fastq) < 2, 'Multiple --fastq files are not supported with gaffe')
    
    require(options.fastq is None or len(options.fastq) in [1, 2], 'Exacty 1 or 2 files must be'
            ' passed with --fastq')
    require(options.interleaved == False or options.fastq is None or len(options.fastq) == 1,
            '--interleaved cannot be used when > 1 fastq given')
    require(sum(map(lambda x : 1 if x else 0, [options.fastq, options.gam_input_reads, options.bam_input_reads])) == 1,
            'reads must be speficied with either --fastq or --gam_input_reads or --bam_input_reads')
    require(options.mapper == 'mpmap' or options.snarls_index is None,
            '--snarls_index can only be used with --mapper mpmap') 
    if options.mapper == 'mpmap':
        require('-S' in context.config.mpmap_opts or '--single-path-mode' in context.config.mpmap_opts,
                '-S must be used with mpmap mapper to produce GAM output')
        require(not options.bam_output,
                '--bam_output not currently supported with mpmap mapper')
    require (not options.bam_output or not options.surject,
             '--bam_output cannot be used in combination with --surject')
    require (not options.id_ranges or not options.surject,
             '--surject not currently supported with --id_ranges')
        

def run_pedigree(job, context, fastq_proband, gam_input_reads_proband, bam_input_reads_proband,
                fastq_maternal, gam_input_reads_maternal, bam_input_reads_maternal,
                fastq_paternal, gam_input_reads_paternal, bam_input_reads_paternal,
                fastq_siblings, gam_input_reads_siblings, bam_input_reads_siblings,
                proband_name, maternal_name, paternal_name, siblings_names,
                interleaved, mapper, indexes,
                reads_file_ids_proband=None, reads_file_ids_maternal=None, reads_file_ids_paternal=None, reads_file_ids_siblings=None,
                reads_chunk_ids_proband=None, reads_chunk_ids_maternal=None, reads_chunk_ids_paternal=None, reads_chunk_ids_siblings=None,
                bam_output=False, surject=False, 
                gbwt_penalty=None, validate=False):
    """
    Split the fastq, then align each chunk.
    
    Exactly one of fastq, gam_input_reads, or bam_input_reads should be
    non-falsey, to indicate what kind of data the file IDs in reads_file_ids or
    reads_chunk_ids correspond to.
    
    Exactly one of reads_file_ids or read_chunks_ids should be specified.
    reads_file_ids holds a list of file IDs of non-chunked input read files,
    which will be chunked if necessary. reads_chunk_ids holds lists of chunk
    IDs for each read file, as produced by run_split_reads_if_needed.
    
    indexes is a dict from index type ('xg', 'gcsa', 'lcp', 'id_ranges',
    'gbwt', 'minimizer', 'distance', 'snarls') to index file ID. Some indexes
    are extra and specifying them will change mapping behavior. Some indexes
    are required for certain values of mapper.
    
    mapper can be 'map', 'mpmap', or 'gaffe'. For 'map' and 'mpmap', the 'gcsa'
    and 'lcp' indexes are required. For 'gaffe', the 'gbwt', 'minimizer' and
    'distance' indexes are required. All the mappers require the 'xg' index.
    
    If bam_output is set, produce BAMs. If surject is set, surject reads down
    to paths. 
    
    If the 'gbwt' index is present and gbwt_penalty is specified, the default
    recombination penalty will be overridden.
    
    returns output gams, one per chromosome, the total mapping time (excluding
    toil-vg overhead such as transferring and splitting files), and output
    BAMs, one per chromosome, if computed.
    """
    
    # Make sure we have exactly one type of input
    assert (bool(fastq_proband) + bool(gam_input_reads_proband) + bool(bam_input_reads_proband) == 1)
    assert (bool(fastq_maternal) + bool(gam_input_reads_maternal) + bool(bam_input_reads_maternal) == 1)
    assert (bool(fastq_paternal) + bool(gam_input_reads_paternal) + bool(bam_input_reads_paternal) == 1)
    assert (bool(fastq_siblings) + bool(gam_input_reads_siblings) + bool(bam_input_reads_siblings) <= 1)
    
    # Make sure we have exactly one kind of file IDs
    assert(bool(reads_file_ids_proband) + bool(reads_chunk_ids_proband) == 1)
    assert(bool(reads_file_ids_maternal) + bool(reads_chunk_ids_maternal) == 1)
    assert(bool(reads_file_ids_paternal) + bool(reads_chunk_ids_paternal) == 1)
    assert(bool(reads_file_ids_siblings) + bool(reads_chunk_ids_siblings) <= 1)
    
    # to encapsulate everything under this job
    proband_mapping_calling_child_jobs = Job()
    maternal_mapping_calling_child_jobs = Job()
    paternal_mapping_calling_child_jobs = Job()
    parental_graph_construction_job = Job()
    proband_sibling_mapping_calling_child_jobs = Job()
    pedigree_joint_calling_job = Job()
    proband_sibling_mapping_calling_child_jobs.addFollowOn(pedigree_joint_calling_job)
    parental_graph_construction_job.addFollowOn(proband_sibling_mapping_calling_child_jobs)
    
    job.addChild(proband_mapping_calling_child_jobs)
    job.addChild(maternal_mapping_calling_child_jobs)
    job.addChild(paternal_mapping_calling_child_jobs)
    job.addFollowOn(parental_graph_construction_job)
    
    proband_first_mapping_job = proband_mapping_calling_child_jobs.addChildJobFn(run_mapping, context, fastq_proband,
                                     gam_input_reads_proband, bam_input_reads_proband,
                                     proband_name,
                                     options.interleaved, options.mapper, importer.resolve(indexes),
                                     reads_file_ids_proband,
                                     bam_output=True, surject=True,
                                     validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    # Start the proband alignment and variant calling
    proband_calling_output = proband_mapping_calling_child_jobs.addFollowOnJobFn(run_pipeline_call, context, options, indexes['xg'], indexes.get('id_ranges'),
                                chr_bam_ids=proband_first_mapping_job.rv(2), baseline_vcf_id, baseline_tbi_id,
                                fasta_id, bed_id, cores=context.config.misc_cores, memory=context.config.misc_mem,
                                disk=context.config.misc_disk).rv()
    
    maternal_mapping_job = maternal_mapping_calling_child_jobs.addChildJobFn(run_mapping, context, fastq_maternal,
                                     gam_input_reads_maternal, bam_input_reads_maternal,
                                     maternal_name,
                                     options.interleaved, options.mapper, importer.resolve(indexes),
                                     reads_file_ids_maternal,
                                     bam_output=True, surject=True,
                                     validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    # Start the maternal alignment and variant calling
    maternal_calling_output = maternal_mapping_calling_child_jobs.addFollowOnJobFn(run_pipeline_call, context, options, indexes['xg'], indexes.get('id_ranges'),
                                chr_bam_ids=maternal_mapping_job.rv(2), baseline_vcf_id, baseline_tbi_id,
                                fasta_id, bed_id, cores=context.config.misc_cores, memory=context.config.misc_mem,
                                disk=context.config.misc_disk).rv()
    
    paternal_mapping_job = paternal_mapping_calling_child_jobs.addChildJobFn(run_mapping, context, fastq_paternal,
                                     gam_input_reads_paternal, bam_input_reads_paternal,
                                     paternal_name,
                                     options.interleaved, options.mapper, importer.resolve(indexes),
                                     reads_file_ids_paternal,
                                     bam_output=True, surject=True,
                                     validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
    # Start the paternal alignment and variant calling
    paternal_calling_output = paternal_mapping_calling_child_jobs.addFollowOnJobFn(run_pipeline_call, context, options, indexes['xg'], indexes.get('id_ranges'),
                                chr_bam_ids=paternal_mapping_job.rv(2), baseline_vcf_id, baseline_tbi_id,
                                fasta_id, bed_id, cores=context.config.misc_cores, memory=context.config.misc_mem,
                                disk=context.config.misc_disk).rv()
    
    
    proband_mapping_output = proband_mapping_root_job.rv()
    maternal_mapping_output = maternal_mapping_root_job.rv()
    paternal_mapping_output = maternal_mapping_root_job.rv()
   
    #TODO: ADD TRIO VARIANT CALLING HERE
    #       -- incorporate gvcf calling job functions here
    #       -- adapt chunk calling infrastructure here
    #TODO: ADD PARENTAL GRAPH CONSTRUCTION HERE
    #TODO: HOOK PARENTAL GRAPH INTO PROBAND AND SIBLING ALIGNMENT HERE
    sibling_root_job_dict =  {}
    if siblings_names is not None: 
        for sibling_number in xrange(len(siblings_names)):
            
            reads_file_ids_siblings_list = []
            if fastq_siblings:
                reads_file_ids_siblings_list = reads_file_ids_siblings[sibling_number*2:(sibling_number*2)+2]
            elif gam_input_reads_siblings or bam_input_reads_siblings:
                reads_file_ids_siblings_list = reads_file_ids_siblings[sibling_number]
            sibling_root_job_dict[sibling_number] = Job.wrapJobFn(run_mapping, context, fastq_siblings[sibling_number],
                                             gam_input_reads_siblings[sibling_number], bam_input_reads_siblings[sibling_number],
                                             siblings_names[sibling_number],
                                             options.interleaved, options.mapper, importer.resolve(indexes),
                                             reads_file_ids_siblings_list,
                                             bam_output=True, surject=True,
                                             validate,
                                             cores=context.config.misc_cores,
                                             memory=context.config.misc_mem,
                                             disk=context.config.misc_disk)
            # Start the sibling alignment
            job.addChild(sibling_root_job_dict[sibling_number])
    
    sibling_mapping_output_dict = {}
    if siblings_names is not None:
        for sibling_number in xrange(len(siblings_names)):
            sibling_mapping_output_dict[sibling_number] = sibling_root_job_dict[sibling_number].rv()
    
    return align_job.rv()


def pedigree_main(context, options):
    """
    Wrapper for vg pedigree. 
    """

    validate_pedigree_options(context, options)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            importer = AsyncImporter(toil)
            
            # Make an index collection
            indexes = {}
           
            # Upload each index we have
            if options.xg_index is not None:
                indexes['xg'] = importer.load(options.xg_index)
            if options.gcsa_index is not None:
                indexes['gcsa'] = importer.load(options.gcsa_index)
                indexes['lcp'] = importer.load(options.gcsa_index + ".lcp")
            if options.gbwt_index is not None:
                indexes['gbwt'] = importer.load(options.gbwt_index)
            if options.distance_index is not None:
                indexes['distance'] = importer.load(options.distance_index)
            if options.minimizer_index is not None:
                indexes['minimizer'] = importer.load(options.minimizer_index)
            if options.snarls_index is not None:
                indexes['snarls'] = importer.load(options.snarls_index)
            if options.id_ranges is not None:
                indexes['id_ranges'] = importer.load(options.id_ranges)
            
            # Upload other local files to the remote IO Store
            inputReadsFileIDsProband = []
            if options.fastq_proband:
                for sample_reads in options.fastq_proband:
                    inputReadsFileIDsProband.append(importer.load(sample_reads))
            elif options.gam_input_reads_proband:
                inputReadsFileIDsProband.append(importer.load(options.gam_input_reads_proband))
            else:
                assert options.bam_input_reads_proband
                inputReadsFileIDsProband.append(importer.load(options.bam_input_reads_proband))

            inputReadsFileIDsMaternal = []
            if options.fastq_maternal:
                for sample_reads in options.fastq_maternal:
                    inputReadsFileIDsMaternal.append(importer.load(sample_reads))
            elif options.gam_input_reads_maternal:
                inputReadsFileIDsMaternal.append(importer.load(options.gam_input_reads_maternal))
            else:
                assert options.bam_input_reads_maternal
                inputReadsFileIDsMaternal.append(importer.load(options.bam_input_reads_maternal))
            
            inputReadsFileIDsPaternal = []
            if options.fastq_paternal:
                for sample_reads in options.fastq_paternal:
                    inputReadsFileIDsPaternal.append(importer.load(sample_reads))
            elif options.gam_input_reads_paternal:
                inputReadsFileIDsPaternal.append(importer.load(options.gam_input_reads_paternal))
            else:
                assert options.bam_input_reads_paternal
                inputReadsFileIDsPaternal.append(importer.load(options.bam_input_reads_paternal))
            
            inputReadsFileIDsSiblings = []
            if options.fastq_siblings:
                for sample_reads in options.fastq_siblings:
                    inputReadsFileIDsSiblings.append(importer.load(sample_reads))
            elif options.gam_input_reads_siblings:
                for sample_gam_reads in options.gam_input_reads_siblings:
                    inputReadsFileIDsSiblings.append(importer.load(sample_gam_reads))
            elif options.bam_input_reads_siblings:
                for sample_bam_reads in options.bam_input_reads_siblings:
                    inputReadsFileIDsSiblings.append(importer.load(sample_bam_reads))
            
            importer.wait()

            # Make a root job
            root_job = Job.wrapJobFn(run_pedigree, context, options.fastq_proband,
                                     options.gam_input_reads_proband, options.bam_input_reads_proband,
                                     options.fastq_maternal,
                                     options.gam_input_reads_maternal, options.bam_input_reads_maternal,
                                     options.fastq_paternal,
                                     options.gam_input_reads_paternal, options.bam_input_reads_paternal,
                                     options.fastq_siblings,
                                     options.gam_input_reads_siblings, options.bam_input_reads_siblings,
                                     options.proband_name,
                                     options.maternal_name,
                                     options.paternal_name,
                                     options.sibling_names,
                                     options.interleaved, options.mapper, importer.resolve(indexes),
                                     reads_file_ids_proband=importer.resolve(inputReadsFileIDsProband),
                                     reads_file_ids_maternal=importer.resolve(inputReadsFileIDsMaternal),
                                     reads_file_ids_paternal=importer.resolve(inputReadsFileIDsPaternal),
                                     reads_file_ids_siblings=importer.resolve(inputReadsFileIDsSiblings),
                                     bam_output=options.bam_output, surject=options.surject,
                                     validate=options.validate,
                                     cores=context.config.misc_cores,
                                     memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)

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
    
