#!/usr/bin/env python2.7
"""
Generate a VCF from a GAM and XG by splitting into GAM/VG chunks.
Chunks are then called in series, and the VCFs stitched together.
Any step whose expected output exists is skipped unles --overwrite 
specified.  

old docker vg tool= quay.io/ucsc_cgl/vg:1.4.0--4cbd3aa6d2c0449730975517fc542775f74910f3
latest docker vg tool=quay.io/ucsc_cgl/vg:latest
"""

import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json
from uuid import uuid4

from toil.job import Job
from toil_lib.toillib import *
from toil_lib.programs import docker_call

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
        
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
    parser.add_argument("--chunk", type=int, default=10000000,
                        help="chunk size")
    parser.add_argument("--overlap", type=int, default=2000,
                        help="amount of overlap between chunks")
    parser.add_argument("--filter_opts", type=str,
                        default="-r 0.9 -afu -s 2 -o 0 -q 15",
                        help="options to pass to vg filter. wrap in \"\"")
    parser.add_argument("--pileup_opts", type=str,
                        default="-w 40 -m 10 -q 10 -a",
                        help="options to pass to vg pileup. wrap in \"\"")
    parser.add_argument("--call_opts", type=str,
                        default="-b 0.4 -f 0.25 -d 10 -s 1",
                        help="options to pass to vg call. wrap in \"\"")
    parser.add_argument("--threads", type=int, default=20,
                        help="number of threads to use in vg call and vg pileup")
    parser.add_argument("--overwrite", action="store_true",
                        help="always overwrite existing files")
                        
    args = args[1:]
        
    return parser.parse_args(args)

def run(cmd, proc_stdout = sys.stdout, proc_stderr = sys.stderr,
        check = True):
    """ run command in shell and throw exception if it doesn't work 
    """
    print cmd
    proc = subprocess.Popen(cmd, shell=True, bufsize=-1,
                            stdout=proc_stdout, stderr=proc_stderr)
    output, errors = proc.communicate()
    sts = proc.wait()
    
    if check is True and sts != 0:
        RealTimeLogger.get().info("Run command error:\nSTDOUT: {}\nSTDERR {}".format(output, errors))
        raise RuntimeError("Command: %s exited with non-zero status %i" % (cmd, sts))
    return output, errors

def merge_call_opts(contig, offset, length, call_opts, sample_name):
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
        elif uo in ["-S", "--sample"]:
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
        opts += " -S {}".format(sample_name)
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

def chunk_gam(gam_path, xg_path, path_name, out_dir, chunks, filter_opts, overwrite):
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
        command = ['filter', os.path.basename(gam_path), '-x', os.path.basename(xg_path), '-R', os.path.basename(chunk_path), '-B', os.path.basename(out_file)] + filter_opts.split(" ")
        docker_call(work_dir=out_dir, parameters=command,
            tool='quay.io/ucsc_cgl/vg:latest')
    
def xg_path_node_id(xg_path, path_name, offset, out_dir):
    """ use vg find to get the node containing a given path position """
    #NOTE: vg find -p range offsets are 0-based inclusive.  
    tmp_out_filename = "{}/tmp_out_{}".format(out_dir, uuid4())
    with open(tmp_out_filename, "w") as tmp_out_file:
        command = ['vg find -x {} -p {}:{}-{}'.format(os.path.basename(xg_path), str(path_name), str(offset), str(offset)), 
                    'vg mod -o -', 'vg view -j -'] 
        docker_call(work_dir=out_dir, parameters=command,
                    tools='quay.io/ucsc_cgl/vg:latest',
                    outfile=tmp_out_file)
        
    command = ['cat data/{}'.format(os.path.basename(tmp_out_filename)), 'jq .node[0].id -']
    stdout = docker_call(work_dir=out_dir, parameters=command,
                    tools='devorbitus/ubuntu-bash-jq-curl',
                    check_output=True)
    
    return int(stdout)

def xg_path_predecessors(xg_path, path_name, node_id, out_dir, context = 1):
    """ get nodes before given node in a path. """
    
    stdout = ''
    command = ['vg find -x {} -n {} -c {}'.format(os.path.basename(xg_path), str(node_id), str(context)),
                'vg view -j -']
    stdout = docker_call(work_dir=out_dir, parameters=command,
                tools='quay.io/ucsc_cgl/vg:latest',
                check_output=True)
    
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

def chunk_vg(xg_path, path_name, out_dir, chunks, chunk_i, overwrite):
    """ use vg find to make one chunk of the graph """
    chunk = chunks[chunk_i]
    vg_chunk_path = chunk_base_name(chunk[0], out_dir, chunk_i, ".vg")
    if overwrite or not os.path.isfile(vg_chunk_path):
        first_node = xg_path_node_id(xg_path, chunk[0], int(chunk[1]), out_dir)
        # xg_path query takes 0-based inclusive coordinates, so we
        # subtract 1 below to convert from BED chunk (0-based exlcusive)
        last_node = xg_path_node_id(xg_path, chunk[0], chunk[2] - 1, out_dir)
        assert first_node > 0 and last_node >= first_node
        # todo: would be cleaner to not have to pad context here
        
        with open(vg_chunk_path, "w") as vg_chunk_path_stream:
            command = ['find', '-x', os.path.basename(xg_path), '-r', str(first_node)+':'+str(last_node), '-c', '1']
            docker_call(work_dir=out_dir, parameters=command,
                        tool='quay.io/ucsc_cgl/vg:latest',
                        outfile=vg_chunk_path_stream)
        
        # but because we got a context, manually go in and make sure
        # our path starts at first_node by deleting everything before
        left_path_padding = xg_path_predecessors(xg_path, path_name, first_node,
                                                 out_dir, context = 1)
        for destroy_id in left_path_padding:
            # destroy should take node list
            destroy_list = vg_chunk_path + ".destroy"

            with open(destroy_list, "w") as destroy_list_stream:
                command = ['vg mod -y {} {}'.format(str(destroy_id), os.path.basename(vg_chunk_path)),
                            'vg mod -o -']
                docker_call(work_dir=out_dir, parameters=command,
                            tools='quay.io/ucsc_cgl/vg:latest',
                            outfile=destroy_list_stream)
            
            run("mv {} {}".format(
                vg_chunk_path + ".destroy", vg_chunk_path))
          
def xg_path_node_offset(xg_path, path_name, offset, out_dir):
    """ get the offset of the node containing the given position of a path
    """
    # first we find the node
    node_id = xg_path_node_id(xg_path, path_name, offset, out_dir)

    # now we find the offset of the beginning of the node
    command = ['find', '-x', os.path.basename(xg_path), '-P', str(path_name), '-n', str(node_id)]
    stdout = docker_call(work_dir=out_dir, parameters=command, 
                    tool='quay.io/ucsc_cgl/vg:latest',
                    check_output=True) 
   
    toks = stdout.split()
    # if len > 2 then we have a cyclic path, which we're assuming we don't
    assert len(toks) == 2
    assert toks[0] == str(node_id)
    node_offset = int(toks[1])
    # node_offset must be before
    RealTimeLogger.get().info("INFO node_offset:{} offset:{}".format(node_offset, offset))
    assert node_offset <= offset
    # sanity check (should really use node size instead of 1000 here)
    assert offset - node_offset < 1000

    return node_offset
    
def sort_vcf(vcf_path, sorted_vcf_path):
    """ from vcflib """
    run("head -10000 {} | grep \"^#\" > {}".format(
        vcf_path, sorted_vcf_path))
    run("grep -v \"^#\" {} | sort -k1,1d -k2,2n >> {}".format(
        vcf_path, sorted_vcf_path))
   
def call_chunk(job, options, index_dir_id, xg_path, path_name, chunks, chunk_i, path_size, overlap,
               pileup_opts, call_options, sample_name, threads, overwrite):
   
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

    # a chunk can be empty if nothing aligns there.
    if not os.path.isfile(gam_path):
        sys.stderr.write("Warning: chunk not found: {}\n".format(gam_path))
        return
    
    # do the pileup.  this is the most resource intensive step,
    # especially in terms of mermory used.
    pu_path = chunk_base_name(path_name, out_dir, chunk_i, ".pu")
    if overwrite or not os.path.isfile(pu_path):
        with open(pu_path, "w") as pu_path_stream:
            command = ['pileup', os.path.basename(vg_path), os.path.basename(gam_path), '-t', str(threads)] + pileup_opts.split(" ")
            docker_call(work_dir=out_dir, parameters=command,
                        tool='quay.io/ucsc_cgl/vg:latest',
                        outfile=pu_path_stream)

    # do the calling.
    # requires the latest version of vg as of 8/2/2016
    vcf_path = chunk_base_name(path_name, out_dir, chunk_i, ".vcf")
    if overwrite or not os.path.isfile(vcf_path + ".gz"):
        offset = xg_path_node_offset(xg_path, chunk[0], chunk[1], out_dir)
        merged_call_opts = merge_call_opts(chunk[0], offset, path_size,
                                           call_options, sample_name)
        with open(vcf_path + ".us", "w") as vgcall_stdout, open(vcf_path + ".call_log", "w") as vgcall_stderr:
            command=['call', os.path.basename(vg_path), os.path.basename(pu_path), '-t',
                     str(threads), str(merged_call_opts)]
            docker_call(work_dir=out_dir, parameters=command,
                        tool='quay.io/ucsc_cgl/vg:latest',
                        outfile=vgcall_stdout, errfile=vgcall_stderr)
 
        sort_vcf(vcf_path + ".us", vcf_path)
        run("rm {}".format(vcf_path + ".us"))
        command=['bgzip', '{}'.format(os.path.basename(vcf_path))]
        docker_call(work_dir=out_dir, parameters=command,
                    tool='quay.io/cmarkello/htslib:latest')
        command=['tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
        docker_call(work_dir=out_dir, parameters=command,
                    tool='quay.io/cmarkello/htslib:latest')

    # do the vcf clip
    left_clip = 0 if chunk_i == 0 else overlap / 2
    right_clip = 0 if chunk_i == len(chunks) - 1 else overlap / 2
    clip_path = chunk_base_name(path_name, out_dir, chunk_i, "_clip.vcf")
    if overwrite or not os.path.isfile(clip_path):
        with open(clip_path, "w") as clip_path_stream:
            command=['bcftools', 'view', '-r', '{}:{}-{}'.format(path_name, chunk[1] + left_clip + 1, chunk[2] - right_clip), '{}'.format(os.path.basename(vcf_path + ".gz"))]
            docker_call(work_dir=out_dir, parameters=command,
                        tool='quay.io/cmarkello/bcftools:latest',
                        outfile=clip_path_stream)

    # Save clip.vcf files to the output store
    out_store.write_output_file(clip_path, os.path.basename(clip_path))
    
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
                    run("cat {} > {}".format(clip_path, vcf_path))
                    first = False
                else:
                    # add on everythin but header
                    run("grep -v \"^#\" {} >> {}".format(clip_path, vcf_path), check=False)

    # add a compressed indexed version
    if overwrite or not os.path.isfile(vcf_path + ".gz"):
        vcf_gz_file = vcf_path + ".gz"
        with open(vcf_gz_file, "w") as vcf_gz_file_stream:
            command=['bgzip', '-c', '{}'.format(os.path.basename(vcf_path))]
            docker_call(work_dir=out_dir, parameters=command,
                        tool='quay.io/cmarkello/htslib:latest',
                        outfile=vcf_gz_file_stream)
        command=['bcftools', 'tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
        docker_call(work_dir=out_dir, parameters=command,
                    tool='quay.io/cmarkello/bcftools:latest')

    # Save merged vcf files to the output store
    out_store.write_output_file(vcf_path, os.path.basename(vcf_path))
    out_store.write_output_file(vcf_path+".gz", os.path.basename(vcf_path+".gz"))
    out_store.write_output_file(vcf_path+".gz.tbi", os.path.basename(vcf_path+".gz.tbi"))
    
    return os.path.basename(vcf_path)

