#!/usr/bin/env python2.7
"""
vg_evaluation_pipeline.py: Run the mapping and variant calling evaluation on all the servers in
parallel.

BAM files with reads must have been already downloaded.

#### DOCKER TEST
#docker run -it --entrypoint /bin/bash --log-driver=none -v /home/cmarkello/debug_eval_output:/data --rm quay.io/ucsc_cgl/vg:latest -c 'vg map -f input.fq -i -M2 -a -u 0 -U -t 6 graph.vg -x graph.vg.xg -g graph.vg.gcsa -n5 > test.out6.gam'

example run: ./vg_evaluation_pipeline.py --realTimeLogging --logError --logDebug --edge_max 5 --kmer_size 16 --index_mode gcsa-mem --include_primary 'azure:hgvm:hgvmeval-jobstore33' '/home/cmarkello/debug_eval_input/BRCA1.vg' 'ref' 81189 '/home/cmarkello/debug_eval_input/BRCA1/NA12877/NA12877.bam.fq' 'NA12877' '/home/cmarkello/debug_eval_output' 'azure:hgvm:hgvmdebugtest-input' 'azure:hgvm:hgvmdebugtest-output'

example run: ./vg_evaluation_pipeline.py --batchSystem mesos --mesosMaster 10.0.0.5:5050 --realTimeLogging --logError --logDebug --edge_max 5 --kmer_size 16 --index_mode gcsa-mem --include_primary 'azure:hgvm:hgvmeval-jobstore33' --path_name 'ref' --path_size 81189 '/home/cmarkello/debug_eval_input/BRCA1.vg' '/home/cmarkello/debug_eval_input/BRCA1/NA12877/NA12877.bam.fq' 'NA12877' '/home/cmarkello/debug_eval_output' 'azure:hgvm:hgvmdebugtest-input' 'azure:hgvm:hgvmdebugtest-output'

example run: ./vg_evaluation_pipeline.py --batchSystem mesos --mesosMaster 10.0.0.5:5050 --realTimeLogging --logError --logDebug --edge_max 5 --kmer_size 16 --index_mode gcsa-mem --include_primary 'azure:hgvm:hgvmeval-jobstore33' '/home/cmarkello/debug_eval_input/BRCA1.vg' '/home/cmarkello/debug_eval_input/BRCA1/NA12877/NA12877.bam.fq' 'NA12877' '/home/cmarkello/debug_eval_output' 'azure:hgvm:hgvmdebugtest-input' 'azure:hgvm:hgvmdebugtest-output'

old docker vg tool=1.4.0--4cbd3aa6d2c0449730975517fc542775f74910f3
new docker vg tool=latest

chr_length_list obtained from Mike Lin's vg dnanexus pipeline configuration
    chr_label_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    chr_length_list = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import ntpath
import getpass
import pdb

from math import ceil
from subprocess import Popen, PIPE
from Bio import SeqIO

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_lib.programs import docker_call
from toil_vg.dockered_chunked_call import *

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
    parser.add_argument("vg_graph", type=str,
        help="Input vg graph file path")
    parser.add_argument("sample_reads", type=str,
        help="Path to sample reads in fastq format")
    parser.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser.add_argument("out_dir", type=str,
        help="directory where all output will be written")
    parser.add_argument("input_store",
        help="sample input IOStore where input files will be temporarily uploaded")
    parser.add_argument("out_store",
        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")
    parser.add_argument("--path_name", nargs='+', type=str,
        help="Name of reference path in the graph (eg. ref or 17)")
    parser.add_argument("--path_size", nargs='+', type=int,
        help="Size of the reference path in the graph")
    parser.add_argument("--offset", type=int,
        help="chromosomal position offset. e.g. 43044293")
    parser.add_argument("--edge_max", type=int, default=5,
        help="maximum edges to cross in index")
    parser.add_argument("--kmer_size", type=int, default=16,
        help="size of kmers to use in indexing and mapping")
    parser.add_argument("--overwrite", default=False, action="store_true",
        help="overwrite existing result files")
    parser.add_argument("--num_fastq_chunks", type=int, default=3,
        help="number of chunks to split the input fastq file records")
    parser.add_argument("--call_chunk_size", type=int, default=10000000,
        help="chunk size for variant calling")
    parser.add_argument("--restat", default=False, action="store_true",
        help="recompute and overwrite existing stats files")
    parser.add_argument("--reindex", default=False, action="store_true",
        help="don't re-use existing indexed graphs")
    parser.add_argument("--index_mode", choices=["rocksdb", "gcsa-kmer",
        "gcsa-mem"], default="gcsa-mem",
        help="type of vg index to use for mapping")
    parser.add_argument("--include_pruned", action="store_true",
        help="use the pruned graph in the index")
    parser.add_argument("--include_primary", action="store_true",
        help="use the primary path in the index")
    parser.add_argument("--overlap", type=int, default=2000,
        help="overlap option that is passed into make_chunks and call_chunk")
    parser.add_argument("--filter_opts", type=str,
        default="-r 0.9 -d 0.05 -e 0.05 -afu -s 1000 -o 10",
        help="options to pass to chunk_gam. wrap in \"\"")
    parser.add_argument("--pileup_opts", type=str,
        default="-w 40 -m 10 -q 10",
        help="options to pass to vg pileup. wrap in \"\"")
    parser.add_argument("--call_opts", type=str,
        default="-b 0.4 -f 0.25 -d 10 -s 1",
        help="options to pass to vg call. wrap in \"\"")


    return parser.parse_args()

# Reverse complement needs a global translation table
reverse_complement_translation_table = string.maketrans("ACGTN", "TGCAN")
def reverse_complement(sequence):
    """  
    Compute the reverse complement of a DNA sequence.
    
    Follows algorithm from <http://stackoverflow.com/a/26615937>
    """
    
    if isinstance(sequence, unicode):
        # Encode the sequence in ASCII for easy translation
        sequence = sequence.encode("ascii", "replace")
    
    # Translate and then reverse
    return sequence.translate(reverse_complement_translation_table)[::-1]
    
def count_Ns(sequence):
    """  
    Return the number of N bases in the given DNA sequence
    """
    
    n_count = 0
    for item in sequence:
        if item == "N": 
            n_count += 1 
     
    return n_count

def get_files_by_file_size(dirname, reverse=False):
    """ Return list of file paths in directory sorted by file size """

    # Get list of files
    filepaths = []
    for basename in os.listdir(dirname):
        filename = os.path.join(dirname, basename)
        if os.path.isfile(filename):
            filepaths.append(filename)

    # Re-populate list with filename, size tuples
    for i in xrange(len(filepaths)):
        filepaths[i] = (filepaths[i], os.path.getsize(filepaths[i]))

    return filepaths

def run(cmd, proc_stdout = sys.stdout, proc_stderr = sys.stderr,
        check = True):
    """ run command in shell and throw exception if it doesn't work 
    """
    RealTimeLogger.get().info(cmd)
    proc = subprocess.Popen(cmd, shell=True, bufsize=-1,
                            stdout=proc_stdout, stderr=proc_stderr)
    output, errors = proc.communicate()
    sts = proc.wait()
    if check is True and sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" % (cmd, sts))
    return output, errors

def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    
    From http://biopython.org/wiki/Split_large_file
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

def run_indexing(job, options):
    """
    For each server listed in the server_list tsv, kick off child jobs to
    align and evaluate it.

    """
    
    RealTimeLogger.get().info("Starting indexing...")
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    graph_file = ntpath.basename(options.vg_graph)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    # Download local input files from the remote storage container
    graph_dir = work_dir
    robust_makedirs(graph_dir)
    
    graph_filename = "{}/graph.vg".format(graph_dir)
    graph_file_remote_path = graph_file
    input_store.read_input_file(graph_file_remote_path, graph_filename)
    
    
    # Now run the indexer.
    RealTimeLogger.get().info("Indexing {}".format(options.vg_graph))
            
    if options.index_mode == "rocksdb":
        # Make the RocksDB index
        command = ['index', '-s', '-k', str(options.kmer_size), '-e', str(options.edge_max), '-t', str(job.cores), os.path.basename(graph_filename), os.path.basename('{}/{}.index'.format(graph_dir, graph_file))]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest')

    elif (options.index_mode == "gcsa-kmer" or
        options.index_mode == "gcsa-mem"):
        # We want a GCSA2/xg index. We have to prune the graph ourselves.
        # See <https://github.com/vgteam/vg/issues/286>.

        # What will we use as our temp combined graph file (containing only
        # the bits of the graph we want to index, used for deduplication)?
        to_index_filename = "{}/to_index.vg".format(
            work_dir)

        # Where will we save the kmers?
        kmers_filename = "{}/index.graph".format(
            work_dir)

        
        with open(to_index_filename, "w") as to_index_file:

            if options.include_pruned:

                RealTimeLogger.get().info("Pruning {} to {}".format(
                    graph_filename, to_index_filename))

                # Prune out hard bits of the graph
                # and complex regions
                # and short disconnected chunks
                command = ['vg mod -p -l {} -t {} -e {} {}'.format(str(options.kmer_size), str(job.cores), str(options.edge_max), os.path.basename(graph_filename)),
                            'vg mod -S -l {} -t {} -'.format(str(options.kmer_size * 2), str(job.cores))]
                docker_call(work_dir=work_dir, parameters=command,
                            tools='quay.io/ucsc_cgl/vg:latest',
                            outfile=to_index_file)

            if options.include_primary:

                # Then append in the primary path. Since we don't knoiw what
                # "it's called, we retain "ref" and all the 19", "6", etc paths
                # "from 1KG.

                RealTimeLogger.get().info(
                    "Adding primary path to {}".format(to_index_filename))
            
                RealTimeLogger.get().info(
                    "to_index_file: {}".format(to_index_file))
                # See
                # https://github.com/vgteam/vg/issues/318#issuecomment-215102199

                # Generate all the paths names we might have for primary paths.
                # It should be "ref" but some graphs don't listen
                ref_names = (["ref", "x", "X", "y", "Y", "m", "M"] +
                    [str(x) for x in xrange(1, 23)])

                ref_options = []
                for name in ref_names:
                    # Put each in a -r option to retain the path
                    ref_options.append("-r")
                    ref_options.append(name)

                # Retain only the specified paths (only one should really exist)
                command = ['mod', '-N'] + ref_options + ['-t', str(job.cores), os.path.basename(graph_filename)]
                docker_call(work_dir=work_dir, parameters=command,
                            tool='quay.io/ucsc_cgl/vg:latest',
                            inputs=[graph_filename],
                            outfile=to_index_file)
                
        time.sleep(1)

        # Now we have the combined to-index graph in one vg file. We'll load
        # it (which deduplicates nodes/edges) and then find kmers.
        RealTimeLogger.get().info("Finding kmers in {} to {}".format(
            to_index_filename, kmers_filename))

        # Make the GCSA2 kmers file
        with open(kmers_filename, "w") as kmers_file:
            command = ['vg view -v {}'.format(os.path.basename(to_index_filename)),
                       'vg kmers -g -B -k {} -H 1000000000 -T 1000000001 -t {} -'.format(str(options.kmer_size), str(job.cores))]
            docker_call(work_dir=work_dir, parameters=command,
                        tools='quay.io/ucsc_cgl/vg:latest',
                        outfile=kmers_file)

        time.sleep(1)

        # Where do we put the GCSA2 index?
        gcsa_filename = graph_filename + ".gcsa"

        RealTimeLogger.get().info("GCSA-indexing {} to {}".format(
                kmers_filename, gcsa_filename))

        # Make the gcsa2 index. Make sure to use 3 doubling steps to work
        # around <https://github.com/vgteam/vg/issues/301>
        command = ['index', '-t', str(job.cores), '-i', 
            os.path.basename(kmers_filename), '-g', os.path.basename(gcsa_filename),
            '-X', '3']
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest')

        # Where do we put the XG index?
        xg_filename = graph_filename + ".xg"

        RealTimeLogger.get().info("XG-indexing {} to {}".format(
                graph_filename, xg_filename))
        
        command = ['index', '-t', str(job.cores), '-x', 
            os.path.basename(xg_filename), os.path.basename(graph_filename),
            '-X', '3']
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest')

    else:
        raise RuntimeError("Invalid indexing mode: " + options.index_mode)

    # Define a file to keep the compressed index in, so we can send it to
    # the output store.
    index_dir_tgz = "{}/index.tar.gz".format(
        job.fileStore.getLocalTempDir())

    # Now save the indexed graph directory to the file store. It can be
    # cleaned up since only our children use it.
    RealTimeLogger.get().info("Compressing index of {}".format(
        graph_filename))
    index_dir_id = write_global_directory(job.fileStore, graph_dir,
        cleanup=True, tee=index_dir_tgz)

    # Save it as output
    RealTimeLogger.get().info("Uploading index of {}".format(
        graph_filename))
    index_key = ntpath.basename(index_dir_tgz)
    out_store.write_output_file(index_dir_tgz, index_key)
    RealTimeLogger.get().info("Index {} uploaded successfully".format(
        index_key))

    #Split fastq files
    return job.addChildJobFn(run_split_fastq, options, index_dir_id, work_dir, cores=3, memory="4G", disk="2G").rv()

def run_split_fastq(job, options, index_dir_id, work_dir):
    
    RealTimeLogger.get().info("Starting fastq split and alignment...")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)


    # We need the sample fastq for alignment
    sample_filename = os.path.basename(options.sample_reads)
    fastq_file = "{}/input.fq".format(work_dir)
    input_store.read_input_file(sample_filename, fastq_file)
    
    record_iter = SeqIO.parse(open(fastq_file),"fastq")
    
    # Find number of records per fastq chunk
    p1 = Popen(['cat', fastq_file], stdout=PIPE)
    p2 = Popen(['wc', '-l'], stdin=p1.stdout, stdout=PIPE)
    p1.stdout.close()
    num_records_total = int(p2.communicate()[0]) / 4.0
    num_records_fastq_chunk = ceil(num_records_total / float(options.num_fastq_chunks))
    
    num_chunks = 0
    for chunk_id, batch in enumerate(batch_iterator(record_iter, num_records_fastq_chunk)):
        num_chunks += 1
        chunk_id = chunk_id + 1
        filename = "{}/group_{}.fq".format(work_dir, chunk_id)
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fastq")
        handle.close()
        
        # Upload the fastq file chunk
        filename_key = os.path.basename(filename)
        out_store.write_output_file(filename, filename_key)
        
        RealTimeLogger.get().info("Wrote {} records to {}".format(count, filename))

        #Run graph alignment on each fastq chunk
        job.addChildJobFn(run_alignment, options, filename_key, chunk_id, index_dir_id, work_dir, cores=3, memory="12G", disk="2G")
    
    return job.addFollowOnJobFn(run_merge_gam, options, num_chunks, index_dir_id, work_dir, cores=3, memory="4G", disk="2G").rv()


def run_alignment(job, options, filename_key, chunk_id, index_dir_id, work_dir):

    RealTimeLogger.get().info("Starting alignment on {} chunk {}".format(options.sample_name, chunk_id))
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    # How long did the alignment take to run, in seconds?
    run_time = None
    
    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)
    
    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)

    # We need the sample fastq for alignment
    sample_filename = ntpath.basename(options.sample_reads)
    fastq_file = "{}/group_{}.fq".format(work_dir, chunk_id)
    out_store.read_input_file(filename_key, fastq_file)
    
    # And a temp file for our aligner output
    output_file = "{}/{}_{}.gam".format(work_dir, options.sample_name, chunk_id)


    # Open the file stream for writing
    with open(output_file, "w") as alignment_file:

        # Start the aligner and have it write to the file

        # Plan out what to run
        vg_parts = ['map', '-f', os.path.basename(fastq_file),
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
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest',
                    outfile=alignment_file)
        
        # Mark when it's done
        end_time = timeit.default_timer()
        run_time = end_time - start_time
 
    time.sleep(1)

    RealTimeLogger.get().info("Aligned {}. Process took {} seconds.".format(output_file, run_time))
    
    
    # Upload the alignment
    alignment_file_key = ntpath.basename(output_file)
    out_store.write_output_file(output_file, alignment_file_key)
    

def run_merge_gam(job, options, num_chunks, index_dir_id, work_dir):
    
    RealTimeLogger.get().info("Starting gam merging...")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)
    
    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)

    # Define a temp file for our merged alignent output
    output_merged_gam = "{}/{}.gam".format(work_dir, options.sample_name)
    
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


    #Run alignment stats
    job.addChildJobFn(run_stats, options, index_dir_id, alignment_file_key, work_dir, cores=2, memory="4G", disk="2G")

    # Run variant calling on .gams by chromosome if no path_name or path_size options are set 
    return_value = []
    vcf_file_key_list = [] 
    if options.path_name and options.path_size:
        #Run variant calling
        for chr_label, chr_length in itertools.izip(options.path_name, options.path_size):
            vcf_file_key = job.addChildJobFn(run_calling, options, index_dir_id, alignment_file_key, work_dir, chr_label, chr_length, cores=2, memory="12G", disk="2G").rv()
            vcf_file_key_list.append(vcf_file_key)
        return_value = job.addFollowOnJobFn(run_merge_vcf, options, index_dir_id, work_dir, vcf_file_key_list, cores=3, memory="8G", disk="2G").rv()
    else:
        raise RuntimeError("Invalid or non existant path_name(s) and/or path_size(s): {}, {}".format(path_name, path_size))

    return return_value

def run_merge_vcf(job, options, index_dir_id, work_dir, vcf_file_key_list):

    RealTimeLogger.get().info("Completed gam merging and gam path variant calling.")
    RealTimeLogger.get().info("Starting vcf merging vcf files: {}".format(" ".join(vcf_file_key_list)))
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)
    
    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)

    ##### DEBUGGING STATEMENTS #####
    filepaths = get_files_by_file_size(work_dir)
    for (filename, filesize) in filepaths:
        RealTimeLogger.get().info("Working directory file list: {}, {}".format(
            filename, filesize))

    ##### END DEBUGGING STATEMENTS ##### 
   
    vcf_merging_file_key_list = [] 
    for vcf_file_key in vcf_file_key_list:
        vcf_file = "{}/{}.gz".format(work_dir, vcf_file_key)
        vcf_file_idx = "{}.tbi".format(vcf_file)
        out_store.read_input_file(vcf_file_key+".gz", vcf_file)
        out_store.read_input_file(vcf_file_key+".gz"+ ".tbi", vcf_file_idx)
        vcf_merging_file_key_list.append(os.path.basename(vcf_file))

    vcf_merged_file_key = "" 
    if len(vcf_merging_file_key_list) > 1:
        # merge vcf files
        vcf_merged_file_key = "{}.vcf.gz".format(options.sample_name)
        command=['bcftools', 'concat', '-O', 'z', '-o', os.path.basename(vcf_merged_file_key), ' '.join(vcf_merging_file_key_list)]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/cmarkello/bcftools')
        command=['bcftools', 'tabix', '-f', '-p', 'vcf', os.path.basename(vcf_merged_file_key)]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/cmarkello/bcftools')
    else:
        vcf_merged_file_key = vcf_merging_file_key_list[0]

    # save variant calling results to the output store
    vcf_file = "{}/{}".format(work_dir, vcf_merged_file_key)
    vcf_file_idx = "{}/{}.tbi".format(work_dir, vcf_merged_file_key)

    out_store.write_output_file(vcf_file, vcf_merged_file_key)
    out_store.write_output_file(vcf_file_idx, vcf_merged_file_key + ".tbi")

    
    #Run downloader to download output IO store files to local output directory.
    vcf_file_id = job.fileStore.writeGlobalFile(vcf_file)
    vcf_file_idx_id = job.fileStore.writeGlobalFile(vcf_file_idx) 
    downloadList = [[vcf_file_id, vcf_merged_file_key], [vcf_file_idx_id, vcf_merged_file_key+".tbi"]]

    return downloadList


def run_calling(job, options, index_dir_id, alignment_file_key, work_dir, path_name, path_size):
    
    RealTimeLogger.get().info("Running variant calling on path {} from alignment file {}".format(path_name, alignment_file_key))
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)
    
    # Download the indexed graph to a directory we can use
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)

    # How long did the alignment take to run, in seconds?
    run_time = None

    # We know what the xg file in there will be named
    xg_file = "{}/graph.vg.xg".format(graph_dir)
    
    # Download the alignment
    alignment_file = "{}/{}.gam".format(work_dir, options.sample_name)
    out_store.read_input_file(alignment_file_key, alignment_file)
    
    variant_call_dir = work_dir

    # run chunked_call
    chunks = dockered_chunked_call(job, options, out_store, work_dir, index_dir_id, job.cores, xg_file, alignment_file, path_name, path_size, options.sample_name, variant_call_dir, options.call_chunk_size, options.overlap, options.filter_opts, options.pileup_opts, options.call_opts, options.overwrite)

    vcf_file_key = job.addFollowOnJobFn(merge_vcf_chunks, options, work_dir, index_dir_id, path_name, path_size, chunks, options.overwrite, cores=3, memory="8G", disk="2G").rv()
 
    RealTimeLogger.get().info("Completed variant calling on path {} from alignment file {}".format(path_name, alignment_file_key))

    return vcf_file_key

def dockered_chunked_call(job, options, out_store, work_dir, index_dir_id, threads, xg_path, gam_path, path_name, path_size, sample_name, out_dir, chunk=10000000,
                          overlap=2000, filter_opts="-r 0.9 -d 0.05 -e 0.05 -afu -s 1000 -o 10",
                          pileup_opts="-w 40 -m 10 -q 10", call_opts="-b 0.4 -f 0.25 -d 10 -s 1",
                          overwrite=True):
    
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    # make things slightly simpler as we split overlap
    # between adjacent chunks
    assert overlap % 2 == 0

    # compute overlapping chunks
    chunks = make_chunks(path_name, path_size, chunk, overlap)

    # split the gam in one go
    chunk_gam(gam_path, xg_path, path_name, out_dir,
              chunks, filter_opts, overwrite)

    # call every chunk in series
    for chunk_i, chunk in enumerate(chunks):
        # make the graph chunk
        chunk_vg(xg_path, path_name, out_dir, chunks, chunk_i, overwrite)
        vg_path = chunk_base_name(path_name, out_dir, chunk_i, ".vg")
        gam_path = chunk_base_name(path_name, out_dir, chunk_i, ".gam")
        
        # upload split gam files
        out_store.write_output_file(vg_path, os.path.basename(vg_path))
        out_store.write_output_file(gam_path, os.path.basename(gam_path))
        
        job.addChildJobFn(call_chunk, options, work_dir, index_dir_id, xg_path, path_name,
                           out_dir, chunks, chunk_i,
                           path_size, overlap,
                           pileup_opts, call_opts,
                           sample_name, threads,
                           overwrite, cores=3, memory="12G", disk="2G")
    return chunks

def run_upload(job, options, uploadList):
    """
    Upload and file in uploadList to the remote IO store specified
    in the input_store option.
    """
    
    RealTimeLogger.get().info("Uploading files to IO store")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)   
    
    for file_key in uploadList:
        file_basename = ntpath.basename(file_key[1])
        RealTimeLogger.get().info("Uploading {} to {} on IO store".format(file_key[0], file_basename))
        fi = job.fileStore.readGlobalFile(file_key[0])
        input_store.write_output_file(fi, file_basename)

    return job.addChildJobFn(run_indexing, options, cores=3, memory="12G", disk="2G").rv()

def run_stats(job, options, index_dir_id, alignment_file_key, work_dir):
    """
    If the stats aren't done, or if they need to be re-done, retrieve the
    alignment file from the output store under alignment_file_key and compute the
    stats file, saving it under stats_file_key.
    
    Uses index_dir_id to get the graph, and thus the reference sequence that
    each read is aligned against, for the purpose of discounting Ns.
    
    Can take a run time to put in the stats.

    Assumes that stats actually do need to be computed, and overwrites any old
    stats.

    TODO: go through the proper file store (and cache) for getting alignment
    data.
    
    """

    RealTimeLogger.get().info("Computing stats for {}".format(options.sample_name))

    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    # Download the indexed graph to a directory we can use
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)

    # How long did the alignment take to run, in seconds?
    run_time = None

    # We know what the vg file in there will be named
    graph_file = "{}/graph.vg".format(graph_dir)

    # Load the node sequences into memory. This holds node sequence string by
    # ID.
    node_sequences = {}

    # Read the alignments in in JSON-line format
    read_graph_filename = "{}/read_graph.json".format(graph_dir)
    with open(read_graph_filename, "w") as read_graph:
        command = ['view', '-j', os.path.basename(graph_file)]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest',
                    inputs=[graph_file],
                    outfile=read_graph)

    with open(read_graph_filename, "r") as read_graph:
        for line in read_graph:
            # Parse the graph chunk JSON
            graph_chunk = json.loads(line)
            for node_dict in graph_chunk.get("node", []):
                # For each node, store its sequence under its id. We want to crash
                # if a node exists for which one or the other isn't defined.
                node_sequences[node_dict["id"]] = node_dict["sequence"]

    # Declare local files for everything
    stats_file = "{}/{}_stats.json".format(job.fileStore.getLocalTempDir(), options.sample_name)
    alignment_file = "{}/{}_alignment.gam".format(job.fileStore.getLocalTempDir(), options.sample_name)
    stats_file_key = ntpath.basename(stats_file)

    # Download the alignment
    out_store.read_input_file(alignment_file_key, alignment_file)

    # Read the alignments in in JSON-line format
    read_alignment_filename = "{}/read_alignment.json".format(graph_dir)
    with open(read_alignment_filename, "w") as read_alignment:
        command = ['view', '-aj', os.path.basename(alignment_file)]
        docker_call(work_dir=work_dir, parameters=command,
                    tool='quay.io/ucsc_cgl/vg:latest',
                    inputs=[alignment_file],
                    outfile=read_alignment)
    
    # Count up the stats
    stats = {
        "total_reads": 0,
        "total_mapped": 0,
        "total_multimapped": 0,
        "mapped_lengths": collections.Counter(),
        "unmapped_lengths": collections.Counter(),
        "aligned_lengths": collections.Counter(),
        "primary_scores": collections.Counter(),
        "primary_mismatches": collections.Counter(),
        "primary_indels": collections.Counter(),
        "primary_substitutions": collections.Counter(),
        "secondary_scores": collections.Counter(),
        "secondary_mismatches": collections.Counter(),
        "secondary_indels": collections.Counter(),
        "secondary_substitutions": collections.Counter(),
        "run_time": run_time
    }

    last_alignment = None

    with open(read_alignment_filename, "r") as read_alignment:
        for line in read_alignment:
            # Parse the alignment JSON
            alignment = json.loads(line)

            # How long is this read?
            length = len(alignment["sequence"])

            if alignment.has_key("score"):
                # This alignment is aligned.
                # Grab its score
                score = alignment["score"]

                # Get the mappings
                mappings = alignment.get("path", {}).get("mapping", [])

                # Calculate the exact match bases
                matches = 0

                # And total up the instances of indels (only counting those where
                # the reference has no Ns, and which aren't leading or trailing soft
                # clips)
                indels = 0

                # And total up the number of substitutions (mismatching/alternative
                # bases in edits with equal lengths where the reference has no Ns).
                substitutions = 0

                # What should the denominator for substitution rate be for this
                # read? How many bases are in the read and aligned?
                aligned_length = 0

                for mapping_number, mapping in enumerate(mappings):
                    # Figure out what the reference sequence for this mapping should
                    # be

                    position = mapping.get("position", {})
                    if position.has_key("node_id"):
                        # We actually are mapped to a reference node
                        ref_sequence = node_sequences[position["node_id"]]
     
                        # Grab the offset
                        offset = position.get("offset", 0)

                        if mapping.get("is_reverse", False):
                            # We start at the offset base on the reverse strand.

                            # Add 1 to make the offset inclusive as an end poiint                        
                            ref_sequence = reverse_complement(
                                ref_sequence[0:offset + 1])
                        else:
                            # Just clip so we start at the specified offset
                            ref_sequence = ref_sequence[offset:]

                    else:
                        # We're aligned against no node, and thus an empty reference
                        # sequence (and thus must be all insertions)
                        ref_sequence = ""

                    # Start at the beginning of the reference sequence for the
                    # mapping.
                    index_in_ref = 0

                    # Pull out the edits
                    edits = mapping.get("edit", [])

                    for edit_number, edit in enumerate(edits):
                        # An edit may be a soft clip if it's either the first edit
                        # in the first mapping, or the last edit in the last
                        # mapping. This flag stores whether that is the case
                        # (although to actually be a soft clip it also has to be an
                        # insertion, and not either a substitution or a perfect
                        # match as spelled by the aligner).
                        may_be_soft_clip = ((edit_number == 0 and
                            mapping_number == 0) or
                            (edit_number == len(edits) - 1 and
                            mapping_number == len(mappings) - 1))

                        # Count up the Ns in the reference sequence for the edit. We
                        # get the part of the reference string that should belong to
                        # this edit.
                        reference_N_count = count_Ns(ref_sequence[
                            index_in_ref:index_in_ref + edit.get("from_length", 0)])

                        if edit.get("to_length", 0) == edit.get("from_length", 0):
                            # Add in the length of this edit if it's actually
                            # aligned (not an indel or softclip)
                            aligned_length += edit.get("to_length", 0)

                        if (not edit.has_key("sequence") and
                            edit.get("to_length", 0) == edit.get("from_length", 0)):
                            # The edit has equal from and to lengths, but no
                            # sequence provided.

                            # We found a perfect match edit. Grab its length
                            matches += edit["from_length"]

                            # We don't care about Ns when evaluating perfect
                            # matches. VG already split out any mismatches into non-
                            # perfect matches, and we ignore the N-matched-to-N
                            # case.

                        if not may_be_soft_clip and (edit.get("to_length", 0) !=
                            edit.get("from_length", 0)):
                            # This edit is an indel and isn't on the very end of a
                            # read.
                            if reference_N_count == 0:
                                # Only count the indel if it's not against an N in
                                # the reference
                                indels += 1

                        if (edit.get("to_length", 0) ==
                            edit.get("from_length", 0) and
                            edit.has_key("sequence")):
                            # The edit has equal from and to lengths, and a provided
                            # sequence. This edit is thus a SNP or MNP. It
                            # represents substitutions.

                            # We take as substituted all the bases except those
                            # opposite reference Ns. Sequence Ns are ignored.
                            substitutions += (edit.get("to_length", 0) -
                                reference_N_count)

                            # Pull those Ns out of the substitution rate denominator
                            # as well.
                            aligned_length -= reference_N_count

                        # We still count query Ns as "aligned" when not in indels

                        # Advance in the reference sequence
                        index_in_ref += edit.get("from_length", 0)

                # Calculate mismatches as what's not perfect matches
                mismatches = length - matches

                if alignment.get("is_secondary", False):
                    # It's a multimapping. We can have max 1 per read, so it's a
                    # multimapped read.

                    if (last_alignment is None or
                        last_alignment.get("name") != alignment.get("name") or
                        last_alignment.get("is_secondary", False)):

                        # This is a secondary alignment without a corresponding primary
                        # alignment (which would have to be right before it given the
                        # way vg dumps buffers
                        raise RuntimeError("{} secondary alignment comes after "
                            "alignment of {} instead of corresponding primary "
                            "alignment\n".format(alignment.get("name"),
                            last_alignment.get("name") if last_alignment is not None
                            else "nothing"))

                    # Log its stats as multimapped
                    stats["total_multimapped"] += 1
                    stats["secondary_scores"][score] += 1
                    stats["secondary_mismatches"][mismatches] += 1
                    stats["secondary_indels"][indels] += 1
                    stats["secondary_substitutions"][substitutions] += 1
                else:
                    # Log its stats as primary. We'll get exactly one of these per
                    # read with any mappings.
                    stats["total_mapped"] += 1
                    stats["primary_scores"][score] += 1
                    stats["primary_mismatches"][mismatches] += 1
                    stats["primary_indels"][indels] += 1
                    stats["primary_substitutions"][substitutions] += 1

                    # Record that a read of this length was mapped
                    stats["mapped_lengths"][length] += 1

                    # And that a read with this many aligned primary bases was found
                    stats["aligned_lengths"][aligned_length] += 1

                    # We won't see an unaligned primary alignment for this read, so
                    # count the read
                    stats["total_reads"] += 1

            elif not alignment.get("is_secondary", False):
                # We have an unmapped primary "alignment"

                # Count the read by its primary alignment
                stats["total_reads"] += 1

                # Record that an unmapped read has this length
                stats["unmapped_lengths"][length] += 1

            # Save the alignment for checking for wayward secondaries
            last_alignment = alignment

    with open(stats_file, "w") as stats_handle:
        # Save the stats as JSON
        json.dump(stats, stats_handle)

    # Now send the stats to the output store where they belong.
    out_store.write_output_file(stats_file, stats_file_key)

def run_download(toil, options, downloadList):
    """
    Download and save each file in downloadList to the local directory specified
    in the out_dir option.
    """
    
    RealTimeLogger.get().info("Download list: {}".format(downloadList))
    
    # Create output directory if it doesn't exist
    try:
        os.makedirs(options.out_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: raise
    
    for outputFileID in downloadList:
        file_basename = ntpath.basename(outputFileID[1])
        RealTimeLogger.get().info("Downloading {} from out_store to {} on the local machine".format(
            outputFileID[0], os.path.join(options.out_dir, file_basename)))
        toil.exportFile(outputFileID[0], 'file://'+os.path.join(options.out_dir, file_basename))
    
    return

def main():
    """
    Parses command line arguments and do the work of the program.
    "args" specifies the program arguments, with args[0] being the executable
    name. The return value should be used as the program's exit code.
    """


    RealTimeLogger.start_master()
    options = parse_args() # This holds the nicely-parsed options object
    
    with Toil(options) as toil:
        if not toil.options.restart:
            
            # Upload local files to the remote IO Store
            inputGraphFileID = toil.importFile('file://'+options.vg_graph)
            inputReadsFileID = toil.importFile('file://'+options.sample_reads)
            uploadList = [[inputGraphFileID, options.vg_graph], [inputReadsFileID, options.sample_reads]]
            
            # Make a root job
            root_job = Job.wrapJobFn(run_upload, options, uploadList, cores=2, memory="5G", disk="2G")
            
            # Run the job and store the returned list of output files to download
            outputFileIDList = toil.start(root_job)
        else:
            outputFileIDList = toil.restart()
        
        # Download output files to the local machine that runs this script
        run_download(toil, options, outputFileIDList)
    
    print("All jobs completed successfully")
    
    RealTimeLogger.stop_master()

if __name__ == "__main__" :
    try:
        main()
    except UserError as e:
        print(e.message, file=sys.stderr)
        sys.exit(1)

