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
    parser = argparse.ArgumentParser(prog='vg_evaluation_pipeline', description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    
    parser.add_argument("sample_reads", type=str,
        help="Path to sample reads in fastq format")
    parser.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser.add_argument("gcsa_index", type=str,
        help="Path to tar and gzipped folder containing .gcsa, .gcsa.lcp, .xg and .graph index files and graph.vg and to_index.vg files. This is equivalent to the output found in the 'run_indexing' toil job function of this pipeline.")    
    parser.add_argument("out_dir", type=str,
        help="directory where all output will be written")
    parser.add_argument("input_store",
        help="sample input IOStore where input files will be temporarily uploaded")
    parser.add_argument("out_store",
        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")
    parser.add_argument("--kmer_size", type=int, default=10,
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--index_mode", choices=["rocksdb", "gcsa-kmer",
        "gcsa-mem"], default="gcsa-mem",
        help="type of vg index to use for mapping")

    
    # Add mapping options
    map_parse_args(parser)

    # Add common docker options
    add_docker_tool_parse_args(parser)

    options = parser.parse_args()

    return parser.parse_args()


def map_parse_args(parser):
    """ centralize indexing parameters here """

    parser.add_argument("--num_fastq_chunks", type=int, default=3,
        help="number of chunks to split the input fastq file records")
    parser.add_argument("--alignment_cores", type=int, default=3,
        help="number of threads during the alignment step")


def run_split_fastq(job, options, index_dir_id):
    
    RealTimeLogger.get().info("Starting fastq split and alignment...")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)


    # We need the sample fastq for alignment
    sample_filename = os.path.basename(options.sample_reads)
    fastq_file = "{}/input.fq".format(work_dir)
    input_store.read_input_file(sample_filename, fastq_file)
    
    
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
    for chunk_id in xrange(options.num_fastq_chunks):
        num_chunks += 1
        chunk_id += 1
        chunk_record_iter = itertools.islice(record_iter, 0, num_records_fastq_chunk)
        chunk_filename = "{}/group_{}.fq".format(work_dir, chunk_id)
        count = SeqIO.write(chunk_record_iter, chunk_filename, "fastq")
        
        # Upload the fastq file chunk
        filename_key = os.path.basename(chunk_filename)
        out_store.write_output_file(chunk_filename, filename_key)
        RealTimeLogger.get().info("Wrote {} records to {}".format(count, chunk_filename))
        
        #Run graph alignment on each fastq chunk
        job.addChildJobFn(run_alignment, options, filename_key, chunk_id, index_dir_id, cores=options.alignment_cores, memory="4G", disk="2G")

    return job.addFollowOnJobFn(run_merge_gam, options, num_chunks, index_dir_id, cores=3, memory="4G", disk="2G").rv()


def run_alignment(job, options, filename_key, chunk_id, index_dir_id):

    RealTimeLogger.get().info("Starting alignment on {} chunk {}".format(options.sample_name, chunk_id))
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_dir = work_dir

    read_global_directory(job.fileStore, index_dir_id, graph_dir)
    
    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)

    # We need the sample fastq for alignment
    sample_filename = os.path.basename(options.sample_reads)
    fastq_file = "{}/group_{}.fq".format(work_dir, chunk_id)
    out_store.read_input_file(filename_key, fastq_file)
    
    # And a temp file for our aligner output
    output_file = "{}/{}_{}.gam".format(work_dir, options.sample_name, chunk_id)


    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:

        # Start the aligner and have it write to the file

        # Plan out what to run
        vg_parts = ['vg', 'map', '-f', os.path.basename(fastq_file),
            '-i', '-M2', '-W', '500', '-u', '0', '-U', '-t', str(job.cores), os.path.basename(graph_file)]

        if options.index_mode == "rocksdb":
            vg_parts += ['-d', os.path.basename(graph_file+".index"), '-n3', '-k',
                str(options.kmer_size)]
        elif options.index_mode == "gcsa-kmer":
            # Use the new default context size in this case
            vg_parts += ['-x', os.path.basename(graph_file+ ".xg"), '-g', os.path.basename(graph_file + ".gcsa"),
                '-n5', '-k', str(options.kmer_size)]
        elif options.index_mode == "gcsa-mem":
            # Don't pass the kmer size, so MEM matching is used
            vg_parts += ['-x', os.path.basename(graph_file+ ".xg"), '-g', os.path.basename(graph_file+ ".gcsa"),
                '-n5']
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
    
    
    # Upload the alignment
    alignment_file_key = os.path.basename(output_file)
    out_store.write_output_file(output_file, alignment_file_key)
    

def run_merge_gam(job, options, num_chunks, index_dir_id):
    
    RealTimeLogger.get().info("Starting gam merging...")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)

    # Define a temp file for our merged alignent output
    output_merged_gam = "{}/{}.gam".format(work_dir, options.sample_name)
    
    # Download the chunked alignments from the outstore to the local work_dir
    gam_chunk_filelist = []
    for i in xrange(num_chunks):
        chunk_id = i + 1
        output_file = "{}/{}_{}.gam".format(work_dir, options.sample_name, chunk_id)
        out_store.read_input_file(os.path.basename(output_file), output_file)
        gam_chunk_filelist.append(output_file)
         
    with open(output_merged_gam, "w") as output_merged_gam_handle:
        for output_chunk_gam in gam_chunk_filelist:
            with open(output_chunk_gam, "r") as output_chunk_gam_handle:
                for line in output_chunk_gam_handle:
                    output_merged_gam_handle.write(line)

    # Upload the merged alignment file
    alignment_file_key = os.path.basename(output_merged_gam)
    out_store.write_output_file(output_merged_gam, alignment_file_key)

    return alignment_file_key

def run_only_mapping(job, options, inputIndexFileID):
    """ run mapping logic by itself.  
    """
    
    input_store = IOStore.get(options.input_store)

    #upload the fastq to the input store, where it's expected to be 
    input_store.write_output_file(options.sample_reads, os.path.basename(options.sample_reads))
    
    alignment_file_key = job.addChildJobFn(run_split_fastq, options, inputIndexFileID, cores=3, memory="4G", disk="2G").rv()

    return alignment_file_key

def fetch_output_gam(options, alignment_file_key):
    """ run_calling leaves the gam output the output store.  copy these over
    to the output directory 
    to do - having both outstore and out_dir seems redundant """

    # Create output directory if it doesn't exist
    try:
        os.makedirs(options.out_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: raise

    RealTimeLogger.get().info("Downloading {} to {}".format(alignment_file_key, options.out_dir))

    out_store = IOStore.get(options.out_store)

    # Read them out of the output store
    out_store.read_input_file(alignment_file_key, os.path.join(options.out_dir, alignment_file_key))

def main():
    """
    Wrapper for vg indexing. 
    """

    RealTimeLogger.start_master()
    options = parse_args() # This holds the nicely-parsed options object

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:
            
            # Upload local files to the remote IO Store
            inputIndexFileID = toil.importFile('file://'+options.gcsa_index)
            
            # Make a root job
            root_job = Job.wrapJobFn(run_only_mapping, options, inputIndexFileID,                                     
                                     cores=2, memory="5G", disk="2G")
            
            # Run the job and store the returned list of output files to download
            alignment_file_key = toil.start(root_job)
        else:
            alignment_file_key = toil.restart()
            
        # copy the indexes out of the output store and into options.out_dir
        fetch_output_gam(options, alignment_file_key)

    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    RealTimeLogger.stop_master()

if __name__ == "__main__" :
    try:
        main()
    except Exception as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)

