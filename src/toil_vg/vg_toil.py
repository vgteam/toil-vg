#!/usr/bin/env python2.7
"""
vg_toil.py: Run the mapping and variant calling on all the servers in
parallel using the vg framework.

old docker vg tool=1.4.0--4cbd3aa6d2c0449730975517fc542775f74910f3
new docker vg tool=latest

chr_length_list obtained from Mike Lin's vg dnanexus pipeline configuration
    chr_label_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    chr_length_list = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import getpass
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
from toil_vg.vg_vcfeval import *
from toil_vg.vg_config import *
from toil_vg.vg_sim import *
from toil_vg.vg_mapeval import *

logger = logging.getLogger(__name__)

def parse_args():
    """
    Takes in the command-line arguments list (args), and returns a nice argparse
    result with fields for all the options.
    
    Borrows heavily from the argparse documentation examples:
    <http://docs.python.org/library/argparse.html>
    """

    # Construct the parser (which is stored in parser)
    # Module docstring lives in __doc__
    # See http://python-forum.com/pythonforum/viewtopic.php?f=3&t=36847
    # And a formatter class so our examples in the docstring look good. Isn't it
    # convenient how we already wrapped it to 80 characters?
    # See http://docs.python.org/library/argparse.html#formatter-class
    parser = argparse.ArgumentParser(prog='toil-vg', description=main.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    subparsers = parser.add_subparsers(dest='command')
    
    # Config subparser
    parser_config = subparsers.add_parser('generate-config',
                                          help='Prints default config file')
    config_subparser(parser_config)
    
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the Toil VG DNA-seq pipeline')
    pipeline_subparser(parser_run)

    # Index subparser
    parser_index = subparsers.add_parser('index', help='Runs only the vg indexing')
    index_subparser(parser_index)

    # Map subparser
    parser_map = subparsers.add_parser('map', help='Runs only the vg mapping')
    map_subparser(parser_map)

    # Call subparser
    parser_call = subparsers.add_parser('call', help='Runs only the vg calling')
    call_subparser(parser_call)

    # vcfeval subparser
    parser_vcfeval = subparsers.add_parser('vcfeval', help='Compare two VCFs wit rtg vcfeval')
    vcfeval_subparser(parser_vcfeval)

    # sim subparser
    parser_sim = subparsers.add_parser('sim', help='Simulate reads and alignments')
    sim_subparser(parser_sim)

    # mapeval subparser
    parser_mapeval = subparsers.add_parser('mapeval', help='Compare mapping results')
    mapeval_subparser(parser_mapeval)

    return parser.parse_args()


def pipeline_subparser(parser_run):
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser_run)
    
    # General options
    parser_run.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser_run.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser_run.add_argument("--xg_index", type=str,
        help="Path to xg index (to use instead of generating new one)")    
    parser_run.add_argument("--gcsa_index", type=str,
        help="Path to GCSA index (to use instead of generating new one)")
    parser_run.add_argument("--id_ranges", type=str,
        help="Path to file with node id ranges for each chromosome in BED format.  If not"
                            " supplied, will be generated from --graphs)")
    parser_run.add_argument("--vcfeval_baseline", type=str,
        help="Path to baseline VCF file for comparison (must be bgzipped and have .tbi)")
    parser_run.add_argument("--vcfeval_fasta", type=str,
        help="Path to DNA sequence file, required for vcfeval. Maybe be gzipped")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser_run)
    
    # Add common indexing options shared with vg_index
    index_parse_args(parser_run)

    # add common mapping options shared with vg_map
    map_parse_args(parser_run)
    
    # Add common calling options shared with vg_call
    chunked_call_parse_args(parser_run)

    # Add common calling options shared with vg_vcfeval
    vcfeval_parse_args(parser_run)

    # Add common docker options
    add_container_tool_parse_args(parser_run)


# Below are the top level jobs of the toil_vg pipeline.  They
# form a chain of "follow-on" jobs as follows:
#
# run_pipeline_upload --> run_pipeline_index --> run_pipeline_map --> run_pipeline_call
#
# Each part of this chain spawns a child job that can also be run independently
# using one of vg_index.main, vg_map.main etc.
#
# Data is communicated across the chain via the output store (at least for now). 


def run_pipeline_index(job, options, inputGraphFileIDs, inputReadsFileIDs, inputXGFileID,
                       inputGCSAFileID, inputLCPFileID, inputIDRangesFileID,
                       inputVCFFileID, inputTBIFileID,
                       inputFastaFileID, inputBeDFileID):
    """
    All indexing.  result is a tarball in thie output store.  Will also do the fastq
    splitting, which doesn't depend on indexing. 
    """

    if inputXGFileID is None:
        xg_file_id = job.addChildJobFn(run_xg_indexing, options, inputGraphFileIDs,
                                       cores=options.xg_index_cores, memory=options.xg_index_mem,
                                       disk=options.xg_index_disk).rv()
    else:
        xg_file_id = inputXGFileID
        
    if inputGCSAFileID is None:
        gcsa_and_lcp_ids = job.addChildJobFn(run_gcsa_prep, options, inputGraphFileIDs,
                                             cores=options.misc_cores, memory=options.misc_mem,
                                             disk=options.misc_disk).rv()
    else:
        assert inputLCPFileID is not None
        gcsa_and_lcp_ids = inputGCSAFileID, inputLCPFileID

    if inputIDRangesFileID is not None:
        id_ranges_file_id = inputIDRangesFileID
    elif len(inputGraphFileIDs) > 1:
        id_ranges_file_id = job.addChildJobFn(run_id_ranges, options, inputGraphFileIDs,
                                               cores=options.misc_cores,
                                               memory=options.misc_mem, disk=options.misc_disk).rv()
    else:
        # don't bother making id ranges if only one input graph
        id_ranges_file_id = None        

    if not options.single_reads_chunk:
        fastq_chunk_ids = job.addChildJobFn(run_split_reads, options, inputReadsFileIDs,
                                            cores=options.misc_cores, memory=options.misc_mem,
                                            disk=options.misc_disk).rv()
    else:
        RealtimeLogger.info("Bypassing reads splitting because --single_reads_chunk enabled")
        fastq_chunk_ids = [inputReadsFileIDs]

    return job.addFollowOnJobFn(run_pipeline_map, options, xg_file_id, gcsa_and_lcp_ids,
                                id_ranges_file_id, fastq_chunk_ids, inputVCFFileID, inputTBIFileID,
                                inputFastaFileID, inputBeDFileID, cores=options.misc_cores,
                                memory=options.misc_mem, disk=options.misc_disk).rv()

def run_pipeline_map(job, options, xg_file_id, gcsa_and_lcp_ids, id_ranges_file_id, fastq_chunk_ids,
                     baseline_vcf_id, baseline_tbi_id, fasta_id, bed_id):
    """ All mapping, then gam merging.  fastq is split in above step"""

    chr_gam_ids = job.addChildJobFn(run_whole_alignment, options,
                                    xg_file_id, gcsa_and_lcp_ids, id_ranges_file_id, fastq_chunk_ids,
                                    cores=options.misc_cores, memory=options.misc_mem,
                                    disk=options.misc_disk).rv()

    return job.addFollowOnJobFn(run_pipeline_call, options, xg_file_id, id_ranges_file_id,
                                chr_gam_ids, baseline_vcf_id, baseline_tbi_id,
                                fasta_id, bed_id, cores=options.misc_cores, memory=options.misc_mem,
                                disk=options.misc_disk).rv()

def run_pipeline_call(job, options, xg_file_id, id_ranges_file_id, chr_gam_ids, baseline_vcf_id,
                      baseline_tbi_id, fasta_id, bed_id):
    """ Run variant calling on the chromosomes in parallel """

    if id_ranges_file_id:
        chroms = [x[0] for x in parse_id_ranges(job, options, id_ranges_file_id)]
    else:
        chroms = options.chroms
    assert len(chr_gam_ids) == len(chroms)

    vcf_tbi_wg_id_pair = job.addChildJobFn(run_all_calling, options, xg_file_id, chr_gam_ids, chroms,
                                           cores=options.misc_cores, memory=options.misc_mem,
                                           disk=options.misc_disk).rv()

    # optionally run vcfeval at the very end.  output will end up in the outstore.
    # f1 score will be returned.
    if baseline_vcf_id is not None:
        return job.addFollowOnJobFn(run_vcfeval, options, vcf_tbi_wg_id_pair,
                                    baseline_vcf_id, baseline_tbi_id, fasta_id, bed_id,
                                    cores=options.vcfeval_cores,
                                    memory=options.vcfeval_mem,
                                    disk=options.vcfeval_disk).rv()
                      

def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Toil vg DNA-seq pipeline
    
    DNA-seq fastqs are split, aligned to an indexed vg reference graph, and variant-called using
    the vg toolset.
    
    General usage:
    1. Type "toil-vg generate-config": Produce an editable config file in the current working directory.
    2. Type "toil-vg run": Given input graphs (one per chromosome) and reads (fastq file), produce a graph index (index can also be input), graph alignment (GAM), VCF variant calls, and (optionally) VCF comparison results.
    3. Type "toil-vg index": Produce an index from input graph(s).
    4. Type "toil-vg map": Produce graph alignment (gam) for each chromosome from input reads and index
    5. Type "toil-vg call": Produce VCF from input index and chromosome gam(s)

    Please read the README.md located in the source directory for more documentation

    Structure of vg DNA-Seq Pipeline (per sample)
                  
                  
                       > 3 ---->    > 5 ---->
                      / ..      |  / ..     |
                     /  ..      v /  ..     v
        0 --> 1 --> 2 --> 3 --> 4 --> 5 --> 6 --> 7
                     \  ..      ^ \  ..     ^
                      \ ..      |  \ ..     | 
                       > 3 ---->    > 5 --->

    0 = Upload local input files to remote input fileStore
    1 = Index input vg reference graph
    2 = Shard Reads
    3 = Align reads to reference graph
    4 = Merge read alignments
    5 = Run vg variant calls on chromosome-chunked .gam alignment
    6 = Merge variant call files
    7 = Download output files from remote output fileStore to local output directory
    ================================================================================
    """

    args = parse_args()
    
    # Write out our config file that's necessary for all other subcommands
    if args.command == 'generate-config':
        config_main(args)
        return

    # Else we merge in our config file args with the command line args for
    # whatever subparser we're in
    options = apply_config_file_args(args)

    # Relative outstore paths can end up who-knows-where.  Make absolute.
    if options.out_store[0] == '.':
        options.out_store = os.path.abspath(options.out_store)

    if args.command == 'vcfeval':
        vcfeval_main(options)
        return

    # Make sure the output store exists and is valid, and through
    # a little record of the command line options and version inside it. 
    init_out_store(options, args.command)

    if args.command == 'run':
        pipeline_main(options)
    elif args.command == 'index':
        index_main(options)
    elif args.command == 'map':
        map_main(options)
    elif args.command == 'call':
        call_main(options)
    elif args.command == 'sim':
        sim_main(options)
    elif args.command == 'mapeval':
        mapeval_main(options)
        
    
def pipeline_main(options):
    
    # make the docker runner
    options.drunner = ContainerRunner(
        container_tool_map = get_container_tool_map(options))

    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'pipeline'

    if options.graphs:
        require(len(options.chroms) == len(options.graphs), '--chroms and --graphs must have'
                ' same number of arguments')
        
    if not options.xg_index or not options.gcsa_index or not options.id_ranges:
        require(options.graphs and options.chroms, '--chroms and --graphs must be specified'
                ' unless --xg_index --gcsa_index and --id_ranges used')

    require(options.fastq is None or len(options.fastq) in [1, 2], 'Exacty 1 or 2 files must be'
            ' passed with --fastq')
    require(options.interleaved == False or options.fastq is None or len(options.fastq) == 1,
            '--interleaved cannot be used when > 1 fastq given')
    require((options.fastq and len(options.fastq)) != (options.gam_input_reads is not None),
            'reads must be speficied with either --fastq or --gam_reads')    

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
            inputGraphFileIDs = []
            if options.graphs:
                for graph in options.graphs:
                    inputGraphFileIDs.append(import_to_store(toil, options, graph))
            inputReadsFileIDs = []
            if options.fastq:
                for sample_reads in options.fastq:
                    inputReadsFileIDs.append(import_to_store(toil, options, sample_reads))
            else:
                inputReadsFileIDs.append(import_to_store(toil, options, options.gam_input_reads))
            if options.gcsa_index:
                inputGCSAFileID = import_to_store(toil, options, options.gcsa_index)
                inputLCPFileID = import_to_store(toil, options, options.gcsa_index + ".lcp")
            else:
                inputGCSAFileID = None
                inputLCPFileID = None
            if options.xg_index:
                inputXGFileID = import_to_store(toil, options, options.xg_index)
            else:
                inputXGFileID = None
            if options.id_ranges:
                inputIDRangesFileID = import_to_store(toil, options, options.id_ranges)
            else:
                inputIDRangesFileID = None
            if options.vcfeval_baseline is not None:
                assert options.vcfeval_baseline.endswith('.vcf.gz')
                assert options.vcfeval_fasta is not None
                inputVCFFileID = import_to_store(toil, options, options.vcfeval_baseline)
                inputTBIFileID = import_to_store(toil, options, options.vcfeval_baseline + '.tbi')
                inputFastaFileID = import_to_store(toil, options, options.vcfeval_fasta)
                inputBedFileID = import_to_store(toil, options, options.vcfeval_bed_regions) \
                                 if options.vcfeval_bed_regions is not None else None
            else:
                inputVCFFileID = None
                inputTBIFileID = None
                inputFastaFileID = None
                inputBedFileID = None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_pipeline_index, options, inputGraphFileIDs,
                                     inputReadsFileIDs, inputXGFileID, inputGCSAFileID,
                                     inputLCPFileID, inputIDRangesFileID,
                                     inputVCFFileID, inputTBIFileID,
                                     inputFastaFileID, inputBedFileID,
                                     cores=options.misc_cores, memory=options.misc_mem,
                                     disk=options.misc_disk)

            # Run the job and store
            toil.start(root_job)
        else:
            toil.restart()

    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline

    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))

if __name__ == "__main__" :
    try:
        main()
    except Exception as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)
