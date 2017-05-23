#!/usr/bin/env python2.7
"""
vg_map.py: map to vg graph producing gam for each chrom

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

logger = logging.getLogger(__name__)

def map_subparser(parser):
    """
    Create a subparser for mapping.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser.add_argument("xg_index", type=str,
        help="Path to xg index")    
    parser.add_argument("gcsa_index", type=str,
        help="Path to GCSA index")
    parser.add_argument("out_store",
        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("--id_ranges", type=str, default=None,
        help="Path to file with node id ranges for each chromosome in BED format.")
    parser.add_argument("--kmer_size", type=int,
        help="size of kmers to use in gcsa-kmer mapping mode")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add mapping options
    map_parse_args(parser)

    # Add common docker options
    add_container_tool_parse_args(parser)


def map_parse_args(parser, stand_alone = False):
    """ centralize indexing parameters here """
    parser.add_argument("--fastq", nargs='+', 
                        help="Input fastq (possibly compressed), two are allowed, one for each mate")
    parser.add_argument("--gam_input_reads", default=None,
                        help="Input reads in GAM format")
    parser.add_argument("--single_reads_chunk", action="store_true", default=False,
                        help="do not split reads into chunks")
    parser.add_argument("--reads_per_chunk", type=int,
                        help="number of reads for each mapping job")
    parser.add_argument("--alignment_cores", type=int,
                        help="number of threads during the alignment step")
    parser.add_argument("--index_mode", choices=["gcsa-kmer", "gcsa-mem"],
                        help="type of vg index to use for mapping")
    parser.add_argument("--interleaved", action="store_true", default=False,
                        help="treat fastq as interleaved read pairs.  overrides map-args")
    parser.add_argument("--map_opts", type=str,
                        help="arguments for vg map (wrapped in \"\")")

def run_mapping(job, options, xg_file_id, gcsa_and_lcp_ids, id_ranges_file_id, reads_file_ids):
    """ split the fastq, then align each chunk """

    if not options.single_reads_chunk:
        reads_chunk_ids = job.addChildJobFn(run_split_reads, options, reads_file_ids,
                                            cores=options.misc_cores, memory=options.misc_mem,
                                            disk=options.misc_disk).rv()
    else:
        RealtimeLogger.info("Bypassing reads splitting because --single_reads_chunk enabled")
        reads_chunk_ids = [reads_file_ids]
    
    return job.addFollowOnJobFn(run_whole_alignment, options, xg_file_id, gcsa_and_lcp_ids,
                                id_ranges_file_id, reads_chunk_ids, cores=options.misc_cores,
                                memory=options.misc_mem, disk=options.misc_disk).rv()    

def run_split_reads(job, options, reads_file_ids):
    """
    split either fastq or gam input reads into chunks.  returns list of chunk file id lists
    (one for each input reads file)
    """

    # this is a list of lists: one list of chunks for each input reads file
    reads_chunk_ids = []
    if options.fastq and len(options.fastq):
        for fastq_i, reads_file_id in enumerate(reads_file_ids):
            reads_chunk_ids.append(job.addChildJobFn(run_split_fastq, options, fastq_i, reads_file_id,
                                                     cores=options.fq_split_cores,
                                                     memory=options.fq_split_mem, disk=options.fq_split_disk).rv())
    else:
        assert len(reads_file_ids) == 1
        reads_chunk_ids.append(job.addChildJobFn(run_split_gam_reads, options, reads_file_ids[0],
                                                  cores=options.fq_split_cores,
                                                  memory=options.fq_split_mem, disk=options.fq_split_disk).rv())

    return reads_chunk_ids


def run_split_fastq(job, options, fastq_i, sample_fastq_id):
    
    RealtimeLogger.info("Starting fastq split")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample fastq for alignment
    fastq_path = os.path.join(work_dir, os.path.basename(options.fastq[fastq_i]))
    fastq_gzipped = os.path.splitext(fastq_path)[1] == '.gz'
    read_from_store(job, options, sample_fastq_id, fastq_path)

    # Split up the fastq into chunks

    # Make sure chunk size even in case paired interleaved
    chunk_size = options.reads_per_chunk
    if chunk_size % 2 != 0:
        chunk_size += 1

    # 4 lines per read
    chunk_lines = chunk_size * 4

    # Note we do this on the command line because Python is too slow
    uc_fastq_name = fastq_path if not fastq_gzipped else 'input_reads.fq'
    if fastq_gzipped:
        cmd = [['gzip', '-d', '-c', os.path.basename(fastq_path)]]
    else:
        cmd = [['cat', os.path.basename(fastq_path)]]

    cmd.append(['split', '-l', str(chunk_lines),
                '--filter=pigz -p {} > $FILE.fq.gz'.format(max(1, int(options.fq_split_cores) - 1)),
                '-', 'fq_chunk.'])

    options.drunner.call(job, cmd, work_dir = work_dir, tool_name='pigz')

    fastq_chunk_ids = []
    for chunk_name in os.listdir(work_dir):
        if chunk_name.endswith('.fq.gz') and chunk_name.startswith('fq_chunk'):
            fastq_chunk_ids.append(write_to_store(job, options, os.path.join(work_dir, chunk_name)))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Split fastq into {} chunks. Process took {} seconds.".format(len(fastq_chunk_ids), run_time))

    return fastq_chunk_ids

def run_split_gam_reads(job, options, gam_reads_file_id):
    """ split up an input reads file in GAM format
    """
    RealtimeLogger.info("Starting gam split")
    start_time = timeit.default_timer()
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample fastq for alignment
    gam_path = os.path.join(work_dir, os.path.basename(options.gam_input_reads))
    read_from_store(job, options, gam_reads_file_id, gam_path)

    # Split up the gam into chunks

    # Make sure chunk size even in case paired interleaved
    chunk_size = options.reads_per_chunk
    if chunk_size % 2 != 0:
        chunk_size += 1

    # 1 line per read
    chunk_lines = chunk_size * 1

    cmd = [['vg', 'view', '-a', os.path.basename(gam_path)]]
    cmd.append(['split', '-l', str(chunk_lines),
                '--filter=vg view -JaG - > $FILE.gam', '-', 'gam_reads_chunk.'])

    options.drunner.call(job, cmd, work_dir = work_dir)

    gam_chunk_ids = []
    for chunk_name in os.listdir(work_dir):
        if chunk_name.endswith('.gam') and chunk_name.startswith('gam_reads_chunk'):
            gam_chunk_ids.append(write_to_store(job, options, os.path.join(work_dir, chunk_name)))
        
    end_time = timeit.default_timer()
    run_time = end_time - start_time
    RealtimeLogger.info("Split gam into {} chunks. Process took {} seconds.".format(len(gam_chunk_ids), run_time))

    return gam_chunk_ids
    
def run_whole_alignment(job, options, xg_file_id, gcsa_and_lcp_ids, id_ranges_file_id, reads_chunk_ids):
    """ align all fastq chunks in parallel
    """
    
    # this will be a list of lists.
    # gam_chunk_file_ids[i][j], will correspond to the jth path (from id_ranges)
    # for the ith gam chunk (generated from fastq shard i)
    gam_chunk_file_ids = []

    for chunk_id, chunk_filename_ids in enumerate(zip(*reads_chunk_ids)):
        #Run graph alignment on each fastq chunk
        gam_chunk_ids = job.addChildJobFn(run_chunk_alignment, options, chunk_filename_ids, chunk_id,
                                          xg_file_id, gcsa_and_lcp_ids, id_ranges_file_id,
                                          cores=options.alignment_cores, memory=options.alignment_mem,
                                          disk=options.alignment_disk).rv()
        gam_chunk_file_ids.append(gam_chunk_ids)

    return job.addFollowOnJobFn(run_merge_gams, options, id_ranges_file_id, gam_chunk_file_ids,
                                cores=options.misc_cores,
                                memory=options.misc_mem, disk=options.misc_disk).rv()


def run_chunk_alignment(job, options, chunk_filename_ids, chunk_id, xg_file_id, gcsa_and_lcp_ids,
                        id_ranges_file_id):

    RealtimeLogger.info("Starting alignment on {} chunk {}".format(options.sample_name, chunk_id))
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    out_store = IOStore.get(options.out_store)

    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_file = os.path.join(work_dir, "graph.vg")

    xg_file = graph_file + ".xg"
    read_from_store(job, options, xg_file_id, xg_file)
    gcsa_file = graph_file + ".gcsa"
    gcsa_file_id = gcsa_and_lcp_ids[0]
    read_from_store(job, options, gcsa_file_id, gcsa_file)
    lcp_file = gcsa_file + ".lcp"
    lcp_file_id = gcsa_and_lcp_ids[1]
    read_from_store(job, options, lcp_file_id, lcp_file)
    
    # We need the sample reads (fastq(s) or gam) for alignment
    reads_files = []
    reads_ext = 'gam' if options.gam_input_reads else 'fq.gz'
    for j, chunk_filename_id in enumerate(chunk_filename_ids):
        reads_file = os.path.join(work_dir, 'reads_chunk_{}_{}.{}'.format(chunk_id, j, reads_ext))
        read_from_store(job, options, chunk_filename_id, reads_file)
        reads_files.append(reads_file)
    
    # And a temp file for our aligner output
    output_file = os.path.join(work_dir, "{}_{}.gam".format(options.sample_name, chunk_id))

    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:

        # Start the aligner and have it write to the file

        # Plan out what to run
        vg_parts = []
        vg_parts += ['vg', 'map',  os.path.basename(graph_file), '-t', str(options.alignment_cores)]
        for reads_file in reads_files:
            input_flag = '-G' if options.gam_input_reads else '-f'
            vg_parts += [input_flag, os.path.basename(reads_file)]
        vg_parts += options.map_opts

        # Override the -i flag in args with the --interleaved command-line flag
        if options.interleaved is True and '-i' not in vg_parts and '--interleaved' not in vg_parts:
            vg_parts += ['-i']
        elif options.interleaved is False and 'i' in vg_parts:
            del vg_parts[vg_parts.index('-i')]
        if options.interleaved is False and '--interleaved' in vg_parts:
            del vg_parts[vg_parts.index('--interleaved')]

        if options.index_mode == "gcsa-kmer":
            # Use the new default context size in this case
            vg_parts += ['-x', os.path.basename(xg_file), '-g', os.path.basename(gcsa_file),
                '-k', str(options.kmer_size)]
        elif options.index_mode == "gcsa-mem":
            # Don't pass the kmer size, so MEM matching is used
            vg_parts += ['-x', os.path.basename(xg_file), '-g', os.path.basename(gcsa_file)]
        else:
            raise RuntimeError("invalid indexing mode: " + options.index_mode)

        RealtimeLogger.info(
            "Running VG for {} against {}: {}".format(options.sample_name, graph_file,
            " ".join(vg_parts)))
        
        # Mark when we start the alignment
        start_time = timeit.default_timer()
        command = vg_parts
        options.drunner.call(job, command, work_dir = work_dir, outfile=alignment_file)
        
        # Mark when it's done
        end_time = timeit.default_timer()
        run_time = end_time - start_time

    RealtimeLogger.info("Aligned {}. Process took {} seconds.".format(output_file, run_time))

    if id_ranges_file_id is not None:
        # Break GAM into multiple chunks at the end. So we need the file
        # defining those chunks.
        id_ranges_file = os.path.join(work_dir, 'id_ranges.tsv')
        read_from_store(job, options, id_ranges_file_id, id_ranges_file)

        # Chunk the gam up by chromosome
        gam_chunks = split_gam_into_chroms(job, work_dir, options, xg_file, id_ranges_file, output_file)

        # Write gam_chunks to store
        gam_chunk_ids = []
        for gam_chunk in gam_chunks:
            gam_chunk_ids.append(write_to_store(job, options, gam_chunk))

        return gam_chunk_ids
        
    else:
        # We can just report one chunk of everything
        return [write_to_store(job, options, output_file)]

def split_gam_into_chroms(job, work_dir, options, xg_file, id_ranges_file, gam_file):
    """
    Create a Rocksdb index then use it to split up the given gam file
    (a local path) into a separate gam for each chromosome.  
    Return a list of filenames.  the ith filename will correspond
    to the ith path in the options.path_name list
    """

    output_index = gam_file + '.index'
    

    index_cmd = ['vg', 'index', '-N', os.path.basename(gam_file),
                 '-d', os.path.basename(output_index), '-t', str(options.gam_index_cores)]
    options.drunner.call(job, index_cmd, work_dir = work_dir)
        
    # Chunk the alignment into chromosomes using the id ranges
    # (note by using id ranges and 0 context and -a we avoid costly subgraph extraction)

    output_bed_name = os.path.join(work_dir, 'output_bed.bed')
    
    # Now run vg chunk on the gam index to get our gams
    # Note: using -a -c 0 -i bypasses subgraph extraction, which is important
    # as it saves a ton of time and memory
    chunk_cmd = ['vg', 'chunk', '-x', os.path.basename(xg_file),
                 '-a', os.path.basename(output_index), '-c', '0',
                 '-r', os.path.basename(id_ranges_file),
                 '-b', os.path.splitext(os.path.basename(gam_file))[0],
                 '-t', str(options.alignment_cores),
                 '-R', os.path.basename(output_bed_name),
                 '-i']
    
    options.drunner.call(job, chunk_cmd, work_dir = work_dir)
    
    # scrape up the vg chunk results into a list of paths to the output gam
    # chunks and return them.  we expect the filenames to be in the 4th column
    # of the output bed from vg chunk. 
    gam_chunks = []
    with open(output_bed_name) as output_bed:
        for line in output_bed:
            toks = line.split('\t')
            if len(toks) > 3:
                gam_chunks.append(os.path.join(work_dir, os.path.basename(toks[3].strip())))
                assert os.path.splitext(gam_chunks[-1])[1] == '.gam'

    return gam_chunks


def run_merge_gams(job, options, id_ranges_file_id, gam_chunk_file_ids):
    """
    Merge together gams, doing each chromosome in parallel
    """
    
    if id_ranges_file_id is not None:
        # Get the real chromosome names
        id_ranges = parse_id_ranges(job, options, id_ranges_file_id)
        chroms = [x[0] for x in id_ranges]
    else:
        # Dump everything in a single default chromosome chunk with a default
        # name
        chroms = ["default"]
    
    chr_gam_ids = []

    for i, chr in enumerate(chroms):
        shard_ids = [gam_chunk_file_ids[j][i] for j in range(len(gam_chunk_file_ids))]
        chr_gam_id = job.addChildJobFn(run_merge_chrom_gam, options, chr, shard_ids,
                                       cores=options.fq_split_cores,
                                       memory=options.fq_split_mem,
                                       disk=options.fq_split_disk).rv()
        chr_gam_ids.append(chr_gam_id)
    
    return chr_gam_ids


def run_merge_chrom_gam(job, options, chr_name, chunk_file_ids):
    """
    Make a chromosome gam by merging up a bunch of gam ids, one 
    for each  shard.  
    """
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    output_file = os.path.join(work_dir, '{}_{}.gam'.format(options.sample_name, chr_name))

    if len(chunk_file_ids) > 1:
        # Would be nice to be able to do this merge with fewer copies.. 
        with open(output_file, 'a') as merge_file:
            for chunk_gam_id in chunk_file_ids:
                tmp_gam_file = os.path.join(work_dir, 'tmp_{}.gam'.format(uuid4()))
                read_from_store(job, options, chunk_gam_id, tmp_gam_file)
                with open(tmp_gam_file) as tmp_f:
                    shutil.copyfileobj(tmp_f, merge_file)
                
        chr_gam_id = write_to_store(job, options, output_file)
    else:
        chr_gam_id = chunk_file_ids[0]
                
    # checkpoint to out store
    if not options.force_outstore or options.tool == 'map':
        if len(chunk_file_ids) == 1:
            read_from_store(job, options, chunk_file_ids[0], output_file)
        write_to_store(job, options, output_file, use_out_store = True)
            
    return chr_gam_id

def map_main(options):
    """
    Wrapper for vg map. 
    """

    # make the docker runner
    options.drunner = ContainerRunner(
        container_tool_map = get_container_tool_map(options))

    require(options.fastq is None or len(options.fastq) in [1, 2], 'Exacty 1 or 2 files must be'
            ' passed with --fastq')
    require(options.interleaved == False or options.fastq is None or len(options.fastq) == 1,
            '--interleaved cannot be used when > 1 fastq given')
    require((options.fastq and len(options.fastq)) != (options.gam_input_reads is not None),
            'reads must be speficied with either --fastq or --gam_reads')
    
    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'map'

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
            inputXGFileID = import_to_store(toil, options, options.xg_index)
            inputGCSAFileID = import_to_store(toil, options, options.gcsa_index)
            inputLCPFileID = import_to_store(toil, options, options.gcsa_index + ".lcp")
            if options.id_ranges is not None:
                inputIDRangesFileID = import_to_store(toil, options, options.id_ranges)
            else:
                inputIDRangesFileID = None
            inputReadsFileIDs = []
            if options.fastq:
                for sample_reads in options.fastq:
                    inputReadsFileIDs.append(import_to_store(toil, options, sample_reads))
            else:
                inputReadsFileIDs.append(import_to_store(toil, options, options.gam_input_reads))
            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_mapping, options,inputXGFileID,
                                     (inputGCSAFileID, inputLCPFileID),
                                     inputIDRangesFileID, inputReadsFileIDs,
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
    
