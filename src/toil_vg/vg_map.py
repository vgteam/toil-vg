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

from math import ceil
from subprocess import Popen, PIPE
from Bio import SeqIO

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_vg.vg_common import *

def map_subparser(parser):
    """
    Create a subparser for mapping.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("sample_reads", type=str,
        help="Path to sample reads in fastq format")
    parser.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser.add_argument("vg_graph", type=str,
        help="Path to vg graph")
    parser.add_argument("xg_index", type=str,
        help="Path to xg index")    
    parser.add_argument("gcsa_index", type=str,
        help="Path to GCSA index")
    parser.add_argument("out_store",
        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")
    parser.add_argument("--kmer_size", type=int,
        help="size of kmers to use in indexing and mapping")
    # these are used for gam merging.  to-do: harmonize these options which are repeated
    # in all the different tools at this point
    parser.add_argument("--path_name", nargs='+', type=str,
        help="Name of reference path in the graph (eg. ref or 17)")
    parser.add_argument("--path_size", nargs='+', type=int,
        help="Size of the reference path in the graph")    

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add mapping options
    map_parse_args(parser)

    # Add common docker options
    add_docker_tool_parse_args(parser)


def map_parse_args(parser, stand_alone = False):
    """ centralize indexing parameters here """

    parser.add_argument("--num_fastq_chunks", type=int,
        help="number of chunks to split the input fastq file records")
    parser.add_argument("--alignment_cores", type=int,
        help="number of threads during the alignment step")
    parser.add_argument("--index_mode", choices=["gcsa-kmer",
        "gcsa-mem"],
        help="type of vg index to use for mapping")        

def run_split_fastq(job, options, graph_file_id, xg_file_id, gcsa_and_lcp_ids, sample_fastq_id):
    
    RealTimeLogger.get().info("Starting fastq split and alignment...")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    out_store = IOStore.get(options.out_store)
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # We need the sample fastq for alignment
    sample_filename = os.path.basename(options.sample_reads)
    fastq_file = "{}/input.fq".format(work_dir)
    read_from_store(job, options, sample_fastq_id, fastq_file)
    
    # Find number of records per fastq chunk
    p1 = Popen(['cat', fastq_file], stdout=PIPE)
    p2 = Popen(['wc', '-l'], stdin=p1.stdout, stdout=PIPE)
    p1.stdout.close()
    num_records_total = int(p2.communicate()[0]) / 4.0
    num_records_fastq_chunk = ceil(num_records_total / float(options.num_fastq_chunks))
 
    # Iterate through records of fastq_file and stream the number of records per fastq
    #   file chunk into a fastq chunk file
  
    num_chunks = 0 
    record_iter = SeqIO.parse(open(fastq_file),"fastq")
    gam_chunk_file_ids = []
    for chunk_id in xrange(options.num_fastq_chunks):
        num_chunks += 1
        chunk_id += 1
        chunk_record_iter = itertools.islice(record_iter, 0, num_records_fastq_chunk)
        chunk_filename = "{}/group_{}.fq".format(work_dir, chunk_id)
        count = SeqIO.write(chunk_record_iter, chunk_filename, "fastq")

        # write the chunk to the jobstore
        chunk_filename_id = write_to_store(job, options, chunk_filename)
        RealTimeLogger.get().info("Wrote {} records to {}".format(count, chunk_filename))
        
        #Run graph alignment on each fastq chunk
        gam_file_id = job.addChildJobFn(run_alignment, options, chunk_filename_id, chunk_id,
                                        graph_file_id, xg_file_id, gcsa_and_lcp_ids,
                                        cores=options.alignment_cores, memory=options.alignment_mem, disk=options.alignment_disk).rv()
        gam_chunk_file_ids.append(gam_file_id)

    return job.addFollowOnJobFn(run_merge_gam, options, gam_chunk_file_ids, xg_file_id,
                                cores=3, memory="4G", disk="2G").rv()


def run_alignment(job, options, chunk_filename_id, chunk_id, graph_file_id, xg_file_id, gcsa_and_lcp_ids):

    RealTimeLogger.get().info("Starting alignment on {} chunk {}".format(options.sample_name, chunk_id))
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    out_store = IOStore.get(options.out_store)

    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_file = os.path.join(work_dir, "graph.vg")
    read_from_store(job, options, graph_file_id, graph_file)

    xg_file = graph_file + ".xg"
    read_from_store(job, options, xg_file_id, xg_file)
    gcsa_file = graph_file + ".gcsa"
    gcsa_file_id = gcsa_and_lcp_ids[0]
    read_from_store(job, options, gcsa_file_id, gcsa_file)
    lcp_file = gcsa_file + ".lcp"
    lcp_file_id = gcsa_and_lcp_ids[1]
    read_from_store(job, options, lcp_file_id, lcp_file)

    # We need the sample fastq for alignment
    fastq_file = os.path.join(work_dir, 'chunk_{}.fq'.format(chunk_id))
    read_from_store(job, options, chunk_filename_id, fastq_file)
    
    # And a temp file for our aligner output
    output_file = os.path.join(work_dir, "{}_{}.gam".format(options.sample_name, chunk_id))

    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:

        # Start the aligner and have it write to the file

        # Plan out what to run
        vg_parts = []
        if hasattr(options, 'vg_map_args'):
            vg_parts += ['vg', 'map', '-f', os.path.basename(fastq_file), os.path.basename(graph_file)]
            vg_parts += options.vg_map_args
        else:
            vg_parts = ['vg', 'map', '-f', os.path.basename(fastq_file),
                        '-i', '-M2', '-W', '500', '-u', '0', '-U',
                        '-O', '-S', '50', '-a', '-t', str(job.cores), os.path.basename(graph_file)]

        if options.index_mode == "gcsa-kmer":
            # Use the new default context size in this case
            vg_parts += ['-x', os.path.basename(xg_file), '-g', os.path.basename(gcsa_file),
                '-n5', '-k', str(options.kmer_size)]
        elif options.index_mode == "gcsa-mem":
            # Don't pass the kmer size, so MEM matching is used
            vg_parts += ['-x', os.path.basename(xg_file), '-g', os.path.basename(gcsa_file), '-n5']
        else:
            raise RuntimeError("invalid indexing mode: " + options.index_mode)

        RealTimeLogger.get().info(
            "Running VG for {} against {}: {}".format(options.sample_name, graph_file,
            " ".join(vg_parts)))
        
        # Mark when we start the alignment
        start_time = timeit.default_timer()
        command = vg_parts
        options.drunner.call(command, work_dir = work_dir, outfile=alignment_file)
        
        # Mark when it's done
        end_time = timeit.default_timer()
        run_time = end_time - start_time

    RealTimeLogger.get().info("Aligned {}. Process took {} seconds.".format(output_file, run_time))
    
    # Send alignment to store
    alignment_file_id = write_to_store(job, options, output_file)
    return alignment_file_id

def run_merge_gam(job, options, chunk_file_ids, xg_file_id):
    
    RealTimeLogger.get().info("Starting gam merging...")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    out_store = IOStore.get(options.out_store)
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download the chunked alignments from the outstore to the local work_dir
    gam_chunk_filelist = []
    for i, chunk_file_id in enumerate(chunk_file_ids):
        chunk_id = i + 1
        output_file = "{}/{}_{}.gam".format(work_dir, options.sample_name, chunk_id)
        read_from_store(job, options, chunk_file_id, output_file)
        gam_chunk_filelist.append(output_file)


    # list of chromsome name, gam file name (no path) pairs
    chr_gam_keys = []

    if options.path_name is not None:
        # Create a bed file of the different chromosomes to pass to vg filter
        # using information passed in via --path_name and --path_size
        bed_file_path = os.path.join(work_dir, "bed_regions.bed")
        with open(bed_file_path, "w") as bed_file:
            assert len(options.path_name) == len(options.path_size)
            for name, size in zip(options.path_name, options.path_size):
                bed_file.write("{}\t0\t{}\n".format(name, size))

        # Download local input files from the remote storage container
        xg_file = os.path.join(work_dir, "graph.vg.xg")
        read_from_store(job, options, xg_file_id, xg_file)

        for i, output_chunk_gam in enumerate(gam_chunk_filelist):        
            # Use vg filter to merge the gam_chunks back together, while dividing them by
            # chromosome.  todo: should we just create the smaller call-chunks right
            # here at the same time? also, can we better parallelize?
            # todo:  clean up options a bit (using mapping cores for filtering) and
            # see about putting lighter weight filters here too (ex identity)
            filter_cmd = ['vg', 'filter', os.path.basename(output_chunk_gam),
                          '-R', os.path.basename(bed_file_path),
                          '-B', options.sample_name, '-t', str(options.alignment_cores),
                          '-x', os.path.basename(xg_file)]
            if i > 0:
                filter_cmd.append('-A')
                
            options.drunner.call(filter_cmd, work_dir = work_dir)


        for chunk_id, name in enumerate(options.path_name):
            # we are relying on convention of vg filter naming output here
            # if that changes, this will break.
            filter_file_path = os.path.join(
                work_dir, "{}-{}.gam".format(options.sample_name, chunk_id))
            assert os.path.isfile(filter_file_path)

            # rename from chunk id to path name
            chr_filterfile_path = os.path.join(
                work_dir, "{}_{}.gam".format(options.sample_name, name))
            os.rename(filter_file_path, chr_filterfile_path)
            chr_gam_keys.append((name, os.path.basename(chr_filterfile_path)))
                
    else:
        # No path information, just append together without splitting into chromosome
        merged_gam_path = os.path.join(work_dir, "{}.gam".format(options.sample_name))
        if len(gam_chunk_filelist) > 0:
            shutil.copy(gam_chunk_filelist[0], merged_gam_path)
            with open(merged_gam_path, 'a') as merge_file:
                for output_chunk_gam in gam_chunk_filelist[1:]:
                    with open(output_chunk_gam) as chunk_file:
                        shutil.copyfileobj(chunk_file, merge_file)                    
        chr_gam_keys = [(None, os.path.basename(merged_gam_path))]

    chr_gam_ids = []
    for name, alignment_key in chr_gam_keys:
        # Upload the merged alignment file to the job store
        gam_file_id = write_to_store(job, options, os.path.join(work_dir, alignment_key))
        # Checkpoint gam to out store
        if not options.force_outstore:
            write_to_store(job, options, os.path.join(work_dir, alignment_key), use_out_store = True)
        chr_gam_ids.append((name, gam_file_id))

    return chr_gam_ids

def map_main(options):
    """
    Wrapper for vg map. 
    """

    RealTimeLogger.start_master()

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))

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
            
            # Upload local files to the remote IO Store
            inputGraphFileID = import_to_store(toil, options, options.vg_graph)
            inputXGFileID = import_to_store(toil, options, options.xg_index)
            inputGCSAFileID = import_to_store(toil, options, options.gcsa_index)
            inputLCPFileID = import_to_store(toil, options, options.gcsa_index + ".lcp")
            sampleFastqFileID = import_to_store(toil, options, options.sample_reads)
            
            # Make a root job
            root_job = Job.wrapJobFn(run_split_fastq, options, inputGraphFileID, inputXGFileID,
                                     (inputGCSAFileID, inputLCPFileID), sampleFastqFileID,
                                     cores=2, memory="5G", disk="2G")
            
            # Run the job and store the returned list of output files to download
            toil.start(root_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    RealTimeLogger.stop_master()

