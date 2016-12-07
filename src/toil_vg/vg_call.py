#!/usr/bin/env python2.7
"""
Generate a VCF from a GAM and XG by splitting into GAM/VG chunks.
Chunks are then called in series, and the VCFs stitched together.
Any step whose expected output exists is skipped unles --overwrite 
specified.  

old docker vg tool= quay.io/ucsc_cgl/vg:1.4.0--4cbd3aa6d2c0449730975517fc542775f74910f3
latest docker vg tool=quay.io/ucsc_cgl/vg:latest
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno
from uuid import uuid4

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_vg.vg_common import *

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("xg_path", type=str,
                        help="input xg file")
    parser.add_argument("gam_path", type=str,
                        help="input alignment")
    parser.add_argument("path_name", type=str,
                        help="name of reference path in graph (ex chr21)")
    parser.add_argument("path_size", type=int,
                        help="size of the reference path in graph")
    parser.add_argument("sample_name", type=str,
                        help="sample name (ex NA12878)")
    parser.add_argument("out_dir", type=str,
                        help="directory where all output will be written")
    parser.add_argument("input_store",
                        help="sample input IOStore where input files will be temporarily uploaded")
    parser.add_argument("out_store",
                        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")    

    # Add common calling options shared with vg_evaluation_pipeline
    chunked_call_parse_args(parser)

    # Add common docker options shared with vg_evaluation pipeline
    add_docker_tool_parse_args(parser)
                        
    return parser.parse_args()

def chunked_call_parse_args(parser):
    """ centralize calling parameters here """
    parser.add_argument("--overlap", type=int, default=2000,
                        help="overlap option that is passed into make_chunks and call_chunk")
    parser.add_argument("--call_chunk_size", type=int, default=10000000,
                        help="chunk size")
    parser.add_argument("--offset", type=int,
                        help="chromosomal position offset. e.g. 43044293")
    parser.add_argument("--filter_opts", type=str,
                        default="-r 0.9 -fu -s 1000 -o 10",
                        help="options to pass to chunk_gam. wrap in \"\"")
    parser.add_argument("--pileup_opts", type=str,
                        default="-w 40 -m 10 -q 10",
                        help="options to pass to vg pileup. wrap in \"\"")
    parser.add_argument("--call_opts", type=str,
                        default="-b 0.4 -f 0.25 -d 10 -s 1",
                        help="options to pass to vg call. wrap in \"\"")
    parser.add_argument("--genotype_opts", type=str,
                        default="",
                        help="options to pass to vg genotype. wrap in \"\"")
    parser.add_argument("--genotype", action="store_true",
                        help="use vg genotype instead of vg call")
    parser.add_argument("--calling_cores", type=int, default=3,
                        help="number of threads during the variant calling step")
    parser.add_argument("--overwrite", action="store_true",
                        help="always overwrite existing files")        
        
def merge_call_opts(contig, offset, length, call_opts, sample_name, sample_flag = '-S'):
    """ combine input vg call  options with generated options, by adding user offset
    and overriding contigs, sample and sequence length"""
    user_opts = call_opts.split()
    user_offset, user_contig, user_ref, user_sample, user_length  = None, None, None, None, None
    for i, uo in enumerate(user_opts):
        if uo in ["-o", "--offset"]:
            user_offset = int(user_opts[i + 1])
            user_opts[i + 1] = str(user_offset + offset)
        elif uo in ["-c", "--contig"]:
            user_contig = user_opts[i + 1]
        elif uo in ["-r", "--ref"]:
            user_ref = user_opts[i + 1]
        elif uo in [sample_flag, "--sample"]:
            user_sample = user_opts[i + 1]
            user_opts[i + 1] = sample_name
        elif uo in ["-l", "--length"]:
            user_length = user_opts[i + 1]
    opts = " ".join(user_opts)
    if user_offset is None:
        opts += " -o {}".format(offset)
    if user_contig is None:
        opts += " -c {}".format(contig)
    if user_ref is None:
        opts += " -r {}".format(contig)  
    if user_sample is None:
        opts += " {} {}".format(sample_flag, sample_name)
    if user_length is None:
        opts += " -l {}".format(length)
    return opts

def make_chunks(path_name, path_size, chunk_size, overlap):
    """ compute chunks as BED (0-based) 3-tuples: ie
    (chr1, 0, 10) is the range from 0-9 inclusive of chr1
    """
    assert chunk_size > overlap
    covered = 0
    chunks = []
    while covered < path_size:
        start = max(0, covered - overlap)
        end = min(path_size, start + chunk_size)
        chunks.append((path_name, start, end))
        covered = end
    return chunks

def chunk_base_name(path_name, out_dir, chunk_i = None, tag= ""):
    """ centralize naming of output chunk-related files """
    bn = os.path.join(out_dir, "{}-chunk".format(path_name))
    if chunk_i is not None:
        bn += "-{}".format(chunk_i)
    return "{}{}".format(bn, tag)

def chunk_gam(drunner, gam_path, xg_path, path_name, out_dir, chunks, filter_opts, overwrite):
    """ use vg filter to chunk up the gam """
    RealTimeLogger.get().info("Starting chunk_gam")
    # make bed chunks
    chunk_path = os.path.join(out_dir, path_name + "_chunks.bed")
    with open(chunk_path, "w") as f:
        for chunk in chunks:
            f.write("{}\t{}\t{}\n".format(chunk[0], chunk[1], chunk[2]))
    # run vg filter on the gam
    stdout = ''
    if overwrite or not any(
            os.path.isfile(chunk_base_name(path_name, out_dir, i, ".gam")) \
               for i in range(len(chunks))):
        
        out_file = os.path.join(out_dir, path_name + "-chunk")
        command = [['vg', 'filter', os.path.basename(gam_path), '-x', os.path.basename(xg_path), '-R', os.path.basename(chunk_path), '-B', os.path.basename(out_file)] + filter_opts.split(" ")]
        drunner.call(command, work_dir=out_dir)
    
def xg_path_node_id(drunner, xg_path, path_name, offset, out_dir):
    """ use vg find to get the node containing a given path position """
    #NOTE: vg find -p range offsets are 0-based inclusive.  
    tmp_out_filename = "{}/tmp_out_{}".format(out_dir, uuid4())
    with open(tmp_out_filename, "w") as tmp_out_file:
        command = [['vg', 'find', '-x', os.path.basename(xg_path), '-p',
                   '{}:{}-{}'.format(str(path_name), str(offset), str(offset))]]
        command.append(['vg', 'mod', '-o', '-'])
        command.append(['vg', 'view', '-j', '-'])
        drunner.call(command, work_dir=out_dir, outfile=tmp_out_file)

    command = [['jq', '.node[0].id', os.path.basename(tmp_out_filename)]]

    # todo : fix this hack:
    # why do we need to do this here? Is it something to do with the jq
    # docker image? all the other tools seem to check /data/ transparently
    if 'jq' in drunner.docker_tool_map:
        command[0][2] = os.path.join('/data', command[0][2])
        
    stdout = drunner.call(command, work_dir=out_dir, check_output=True, outfile=None)
    
    return int(stdout)

def xg_path_predecessors(drunner, xg_path, path_name, node_id, out_dir, context = 1):
    """ get nodes before given node in a path. """
    
    stdout = ''
    command = [['vg', 'find', '-x', os.path.basename(xg_path), '-n', str(node_id), '-c', str(context)]]
    command.append(['vg', 'view', '-j', '-'])
    stdout = drunner.call(command, work_dir=out_dir, check_output=True)
    
    # get our json graph
    j = json.loads(stdout)
    paths = j["path"]
    path = [x for x in paths if x["name"] == path_name][0]
    mappings = path["mapping"]
    assert len(mappings) > 0
    # check that we have a node_mapping
    assert len([x for x in mappings if x["position"]["node_id"] == node_id]) == 1
    # collect mappings that come before
    out_ids = []
    for mapping in mappings:
        if mapping["position"]["node_id"] == node_id:
            break
        out_ids.append(mapping["position"]["node_id"])
    return out_ids

def chunk_vg(drunner, xg_path, path_name, out_dir, chunks, chunk_i, overwrite):
    """ use vg find to make one chunk of the graph """
    chunk = chunks[chunk_i]
    vg_chunk_path = chunk_base_name(chunk[0], out_dir, chunk_i, ".vg")
    if overwrite or not os.path.isfile(vg_chunk_path):
        first_node = xg_path_node_id(drunner, xg_path, chunk[0], int(chunk[1]), out_dir)
        # xg_path query takes 0-based inclusive coordinates, so we
        # subtract 1 below to convert from BED chunk (0-based exlcusive)
        last_node = xg_path_node_id(drunner, xg_path, chunk[0], chunk[2] - 1, out_dir)
        assert first_node > 0 and last_node >= first_node
        # todo: would be cleaner to not have to pad context here
        
        with open(vg_chunk_path, "w") as vg_chunk_path_stream:
            command = [['vg', 'find', '-x', os.path.basename(xg_path), '-r', str(first_node)+':'+str(last_node), '-c', '1']]
            drunner.call(command, work_dir=out_dir, outfile=vg_chunk_path_stream)
            RealTimeLogger.get().info("Output vg chunk size {}".format(
                os.path.getsize(vg_chunk_path)))
        
        # but because we got a context, manually go in and make sure
        # our path starts at first_node by deleting everything before
        left_path_padding = xg_path_predecessors(drunner, xg_path, path_name, first_node,
                                                 out_dir, context = 1)
        for destroy_id in left_path_padding:
            # destroy should take node list
            destroy_list = vg_chunk_path + ".destroy"

            with open(destroy_list, "w") as destroy_list_stream:
                command = [['vg', 'mod', '-y', str(destroy_id), os.path.basename(vg_chunk_path)]]
                command.append(['vg', 'mod', '-o', '-'])
                drunner.call(command, work_dir=out_dir, outfile=destroy_list_stream)
            
            drunner.call(['mv', vg_chunk_path + ".destroy", vg_chunk_path])
          
def xg_path_node_offset(drunner, xg_path, path_name, offset, out_dir):
    """ get the offset of the node containing the given position of a path
    """
    # first we find the node
    node_id = xg_path_node_id(drunner, xg_path, path_name, offset, out_dir)

    # now we find the offset of the beginning of the node
    command = [['vg', 'find', '-x', os.path.basename(xg_path), '-P', str(path_name), '-n', str(node_id)]]
    stdout = drunner.call(command, work_dir=out_dir, check_output=True) 
   
    toks = stdout.split()
    # if len > 2 then we have a cyclic path, which we're assuming we don't
    assert len(toks) == 2
    assert toks[0] == str(node_id)
    node_offset = int(toks[1])
    # node_offset must be before
    assert node_offset <= offset
    # sanity check (should really use node size instead of 1000 here)
    assert offset - node_offset < 1000

    return node_offset
    
def sort_vcf(drunner, vcf_path, sorted_vcf_path):
    """ from vcflib """
    vcf_dir, vcf_name = os.path.split(vcf_path)
    with open(sorted_vcf_path, "w") as outfile:
        drunner.call([['bcftools', 'view', '-h', vcf_name]], outfile=outfile,
                     work_dir=vcf_dir)
    with open(sorted_vcf_path, "a") as outfile:
        drunner.call([['bcftools', 'view', '-H', vcf_name],
                      ['sort', '-k1,1d', '-k2,2n']], outfile=outfile,
                     work_dir=vcf_dir)

def run_vg_call(options, xg_path, vg_path, gam_path, path_name, chunk, chunk_i, path_size, 
                pileup_opts, call_options, sample_name, work_dir, out_dir, vcf_path,
                threads, overwrite):
    """ Create a VCF with vg call """
                        
    # do the pileup.  this is the most resource intensive step,
    # especially in terms of mermory used.
    pu_path = chunk_base_name(path_name, out_dir, chunk_i, ".pu")
    if overwrite or not os.path.isfile(pu_path):
        with open(pu_path, "w") as pu_path_stream:
            command = [['vg', 'pileup', os.path.basename(vg_path), os.path.basename(gam_path), '-t', str(threads)] + pileup_opts.split(" ")]
            options.drunner.call(command, work_dir=out_dir, outfile=pu_path_stream)

    # do the calling.
    # requires the latest version of vg as of 8/2/2016
    if overwrite or not os.path.isfile(vcf_path + ".gz"):
        offset = xg_path_node_offset(options.drunner, xg_path, chunk[0], chunk[1], out_dir)
        merged_call_opts = merge_call_opts(chunk[0], offset, path_size,
                                           call_options, sample_name)
        with open(vcf_path + ".us", "w") as vgcall_stdout, open(vcf_path + ".call_log", "w") as vgcall_stderr:
            command=[['vg', 'call', os.path.basename(vg_path), os.path.basename(pu_path), '-t',
                     str(threads)] + str(merged_call_opts).split()]
            options.drunner.call(command, work_dir=out_dir,
                        outfile=vgcall_stdout, errfile=vgcall_stderr)
 
                
def run_vg_genotype(options, xg_path, vg_path, gam_path, path_name, chunk, chunk_i, path_size, 
                    genotype_options, sample_name, work_dir, out_dir, vcf_path,
                    threads, overwrite):
    """ Create a VCF with vg genotype """

    # Make a gam index
    gam_index_path = os.path.join(out_dir, os.path.basename(gam_path) + ".index")
    if overwrite or not os.path.isfile(gam_index_path):
        command = ['vg', 'index', '-N', os.path.basename(gam_path), '-d', os.path.basename(gam_index_path)]
        options.drunner.call(command, work_dir=work_dir)

    # Do the genotyping
    if overwrite or not os.path.isfile(vcf_path + ".gz"):
        offset = xg_path_node_offset(options.drunner, xg_path, chunk[0], chunk[1], out_dir)
        merged_genotype_opts = merge_call_opts(chunk[0], offset, path_size,
                                           genotype_options, sample_name, sample_flag = '-s')
        with open(vcf_path + ".us", "w") as vgcall_stdout, open(vcf_path + ".call_log", "w") as vgcall_stderr:

            command=[['vg', 'genotype', os.path.basename(vg_path), os.path.basename(gam_index_path),
                      '-t', str(threads), '-v'] + str(merged_genotype_opts).split()]
            options.drunner.call(command, work_dir=out_dir,
                        outfile=vgcall_stdout, errfile=vgcall_stderr)

def call_chunk(job, options, index_dir_id, xg_path, path_name, chunks, chunk_i, path_size, overlap,
               pileup_opts, call_options, genotype_options, sample_name, threads, overwrite):
   
    RealTimeLogger.get().info("Running call_chunk on path {} and chunk {}".format(path_name, chunk_i))
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    RealTimeLogger.get().info("Attempting to set up the IO stores.")
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download local input files from the remote storage container
    graph_dir = work_dir
    read_global_directory(job.fileStore, index_dir_id, graph_dir)

    """ create VCF from a given chunk """
    out_dir = work_dir
    chunk = chunks[chunk_i]
    path_name = chunk[0]
    vg_path = chunk_base_name(path_name, out_dir, chunk_i, ".vg")
    gam_path = chunk_base_name(path_name, out_dir, chunk_i, ".gam")

    # Download the chunked vg and gam files
    RealTimeLogger.get().info("Attempting to read vg_path and gam_path files.")
    out_store.read_input_file(os.path.basename(vg_path), vg_path) 
    out_store.read_input_file(os.path.basename(gam_path), gam_path)

    vcf_path = chunk_base_name(path_name, out_dir, chunk_i, ".vcf")
    
    # Run vg call
    if options.genotype:
        run_vg_genotype(options, xg_path, vg_path, gam_path, path_name, chunk, chunk_i, path_size, 
                        genotype_options, sample_name, work_dir, out_dir,
                        vcf_path, threads, overwrite)
    else:
        run_vg_call(options, xg_path, vg_path, gam_path, path_name, chunk, chunk_i, path_size, 
                    pileup_opts, call_options, sample_name, work_dir, out_dir,
                    vcf_path, threads, overwrite)

    # Sort the output
    sort_vcf(options.drunner, vcf_path + ".us", vcf_path)
    options.drunner.call(['rm', vcf_path + '.us'])
    command=['bgzip', '{}'.format(os.path.basename(vcf_path))]
    options.drunner.call(command, work_dir=out_dir)
    command=['tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
    options.drunner.call(command, work_dir=out_dir)
    
    # do the vcf clip
    left_clip = 0 if chunk_i == 0 else overlap / 2
    right_clip = 0 if chunk_i == len(chunks) - 1 else overlap / 2
    clip_path = chunk_base_name(path_name, out_dir, chunk_i, "_clip.vcf")
    if overwrite or not os.path.isfile(clip_path):
        with open(clip_path, "w") as clip_path_stream:
            command=['bcftools', 'view', '-r', '{}:{}-{}'.format(path_name, chunk[1] + left_clip + 1, chunk[2] - right_clip), '{}'.format(os.path.basename(vcf_path + ".gz"))]
            options.drunner.call(command, work_dir=out_dir, outfile=clip_path_stream)

    # Save clip.vcf files to the output store
    out_store.write_output_file(clip_path, os.path.basename(clip_path))

def run_calling(job, options, index_dir_id, alignment_file_key, path_name, path_size):
    
    RealTimeLogger.get().info("Running variant calling on path {} from alignment file {}".format(path_name, alignment_file_key))
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
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
    chunks = dockered_chunked_call(job, options, out_store, work_dir, index_dir_id, job.cores, xg_file,
                                   alignment_file, path_name, path_size, options.sample_name, variant_call_dir,
                                   options.call_chunk_size, options.overlap, options.filter_opts, options.pileup_opts,
                                   options.call_opts, options.genotype_opts, options.overwrite)

    vcf_file_key = job.addFollowOnJobFn(merge_vcf_chunks, options, index_dir_id, path_name, path_size, chunks, options.overwrite, cores=2, memory="4G", disk="2G").rv()
 
    RealTimeLogger.get().info("Completed variant calling on path {} from alignment file {}".format(path_name, alignment_file_key))

    return vcf_file_key

def dockered_chunked_call(job, options, out_store, work_dir, index_dir_id, threads,
                          xg_path, gam_path, path_name, path_size, sample_name, out_dir, chunk,
                          overlap, filter_opts, pileup_opts, call_opts, genotype_opts,
                          overwrite):
    """
    dockered_chunked_call IS NOT A TOIL JOB FUNCTION. It is just a helper function.
    """
    
    #if not os.path.isdir(out_dir):
    #   os.makedirs(out_dir)

    # make things slightly simpler as we split overlap
    # between adjacent chunks
    assert overlap % 2 == 0

    # compute overlapping chunks
    chunks = make_chunks(path_name, path_size, chunk, overlap)
    
    # split the gam in one go
    chunk_gam(options.drunner, gam_path, xg_path, path_name, work_dir,
              chunks, filter_opts, overwrite)

    # call every chunk in series
    for chunk_i, chunk in enumerate(chunks):
        # make the graph chunk
        chunk_vg(options.drunner, xg_path, path_name, work_dir, chunks, chunk_i, overwrite)
        
        vg_path = chunk_base_name(path_name, work_dir, chunk_i, ".vg")
        gam_path = chunk_base_name(path_name, work_dir, chunk_i, ".gam")
        
        # upload split gam files
        out_store.write_output_file(vg_path, os.path.basename(vg_path))
        out_store.write_output_file(gam_path, os.path.basename(gam_path))
        
        job.addChildJobFn(call_chunk, options, index_dir_id, xg_path, path_name, chunks, chunk_i,
                           path_size, overlap,
                           pileup_opts, call_opts, genotype_opts,
                           sample_name, threads,
                           overwrite, cores="{}".format(threads), memory="4G", disk="2G")
    return chunks

def merge_vcf_chunks(job, options, index_dir_id, path_name, path_size, chunks, overwrite):
    """ merge a bunch of clipped vcfs created above, taking care to 
    fix up the headers.  everything expected to be sorted already """
    
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    input_store = IOStore.get(options.input_store)
    out_store = IOStore.get(options.out_store)   

    # Define work directory for docker calls
    out_dir = job.fileStore.getLocalTempDir()
    
    # Download local input files from the remote storage container
    read_global_directory(job.fileStore, index_dir_id, out_dir)

    vcf_path = os.path.join(out_dir, path_name + ".vcf")
    
    if overwrite or not os.path.isfile(vcf_path):
        first = True
        for chunk_i, chunk in enumerate(chunks):
            clip_path = chunk_base_name(path_name, out_dir, chunk_i, "_clip.vcf")
            # Download clip.vcf file
            out_store.read_input_file(os.path.basename(clip_path), clip_path)

            if os.path.isfile(clip_path):
                if first is True:
                    # copy everything including the header
                    with open(vcf_path, "w") as outfile:
                        options.drunner.call(['cat', os.path.basename(clip_path)], outfile=outfile,
                                             work_dir=out_dir)
                    first = False
                else:
                    # add on everythin but header
                    with open(vcf_path, "a") as outfile:
                        options.drunner.call(['bcftools', 'view', '-H', os.path.basename(clip_path)],
                                             outfile=outfile, work_dir=out_dir)

    # add a compressed indexed version
    if overwrite or not os.path.isfile(vcf_path + ".gz"):
        vcf_gz_file = vcf_path + ".gz"
        with open(vcf_gz_file, "w") as vcf_gz_file_stream:
            command=['bgzip', '-c', '{}'.format(os.path.basename(vcf_path))]
            options.drunner.call(command, work_dir=out_dir, outfile=vcf_gz_file_stream)
        command=['bcftools', 'tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
        options.drunner.call(command, work_dir=out_dir)

    # Save merged vcf files to the output store
    out_store.write_output_file(vcf_path, os.path.basename(vcf_path))
    out_store.write_output_file(vcf_path+".gz", os.path.basename(vcf_path+".gz"))
    out_store.write_output_file(vcf_path+".gz.tbi", os.path.basename(vcf_path+".gz.tbi"))
    
    return os.path.basename(vcf_path)

def run_only_chunked_call(job, options, inputXGFileID, inputGamFileID):
    """ run chunked_call logic by itself.  todo-modularize existing pipeline enough
    so that we don't need to duplicate this code """

    # run_calling wants the index tarred up insaide toilib's
    # read/write_global_directory interface.  so we handle that here
    # todo: fix up run_calling interface so we don't need to do this
    # note: the file name important here as it's what run_calling expects
    indexDir = job.fileStore.getLocalTempDir()
    indexPath = os.path.join(indexDir, "graph.vg.xg")
    job.fileStore.readGlobalFile(inputXGFileID, userPath=indexPath)
    RealTimeLogger.get().info("Compressing input index to {}".format(indexPath + ".gz"))
    indexDirId = write_global_directory(job.fileStore, indexDir,
                                          cleanup=True, tee=indexPath + ".gz")

    inputGamFileID = options.gam_path
    vcf_file_key = job.addChildJobFn(run_calling, options, indexDirId, inputGamFileID,
                                     options.path_name, options.path_size, cores=options.calling_cores,
                                     memory="4G", disk="2G").rv()

    return job.addFollowOnJobFn(wait_for_vcf, vcf_file_key).rv()

def wait_for_vcf(job, vcf_file_key):
    return vcf_file_key

def fetch_output_vcf(options, vcf_file_key):
    """ run_calling leaves the vcf output the output store.  copy these over
    to the output directory 
    to do - having both outstore and out_dir seems redundant """

    # Create output directory if it doesn't exist
    try:
        os.makedirs(options.out_dir)
    except OSError as exception:
        if exception.errno != errno.EEXIST: raise

    RealTimeLogger.get().info("Downloading {} to {}".format(vcf_file_key, options.out_dir))

    # Derive run_calling output files from its return value
    vcf_file = "{}.gz".format(vcf_file_key)
    vcf_file_idx = "{}.gz.tbi".format(vcf_file_key)
    out_store = IOStore.get(options.out_store)

    # Read them out of the output store
    out_store.read_input_file(vcf_file, os.path.join(options.out_dir, vcf_file))
    out_store.read_input_file(vcf_file, os.path.join(options.out_dir, vcf_file_idx))
    
def main():
    """ no harm in preserving command line access to chunked_call for debugging """
    
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
            inputXGFileID = toil.importFile('file://'+options.xg_path)
            inputGamFileID = toil.importFile('file://'+options.gam_path)

            # Make a root job
            root_job = Job.wrapJobFn(run_only_chunked_call, options, inputXGFileID,
                                     inputGamFileID, cores=2, memory="5G", disk="2G")
            
            # Run the job and store the returned list of output files to download
            vcf_file_key = toil.start(root_job)
        else:
            vcf_file_key = toil.restart()

        # copy the vcf out of the output store and into options.out_dir
        fetch_output_vcf(options, vcf_file_key)
                
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
        
        
