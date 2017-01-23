#!/usr/bin/env python2.7
"""
Generate a VCF from a GAM and XG by splitting into GAM/VG chunks.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, timeit, errno, copy
from uuid import uuid4

from toil.common import Toil
from toil.job import Job
from toil_lib.toillib import *
from toil_vg.vg_common import *
from toil_lib import require


def call_subparser(parser):
    """
    Create a subparser for calling.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("xg_path", type=str,
                        help="input xg file")
    parser.add_argument("gam_path", type=str,
                        help="input alignment")
    parser.add_argument("sample_name", type=str,
                        help="sample name (ex NA12878)")
    parser.add_argument("out_store",
                        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")    

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common calling options shared with vg_evaluation_pipeline
    chunked_call_parse_args(parser)

    # Add common docker options shared with vg_evaluation pipeline
    add_docker_tool_parse_args(parser)
                        

def chunked_call_parse_args(parser):
    """ centralize calling parameters here """
    parser.add_argument("--overlap", type=int,
                        help="overlap option that is passed into make_chunks and call_chunk")
    parser.add_argument("--call_chunk_size", type=int,
                        help="chunk size")
    parser.add_argument("--call_opts", type=str,
                        help="arguments to pass to vg call (wrapped in \"\")")
    parser.add_argument("--genotype", action="store_true",
                        help="use vg genotype instead of vg call")
    parser.add_argument("--genotype_opts", type=str,
                        help="arguments to pass to vg genotype (wrapped in \"\")")
    parser.add_argument("--filter_opts", type=str,
                        help="argument to pass to vg filter (wrapped in \"\")")
    parser.add_argument("--calling_cores", type=int,
                        help="number of threads during the variant calling step")
        
def merge_call_opts(contig, offset, length, call_opts, sample_name, sample_flag = '-S'):
    """ combine input vg call  options with generated options, by adding user offset
    and overriding contigs, sample and sequence length"""
    user_opts = copy.deepcopy(call_opts)
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

def run_vg_call(options, xg_path, vg_path, gam_path, path_name, chunk_offset, path_size, work_dir, vcf_path):
    """ Create a VCF with vg call """
                        
    # do the pileup.  this is the most resource intensive step,
    # especially in terms of mermory used.
    # can stream pileup directly to caller. would save disk but worry about memory. 
    pu_path = os.path.join(work_dir, 'chunk_{}_{}.pu'.format(path_name, chunk_offset))
    with open(pu_path, "w") as pu_path_stream:
        command = [['vg', 'filter', os.path.basename(gam_path),
                    '-x', os.path.basename(xg_path), '-t', '1'] + options.filter_opts]
        command.append(['vg', 'pileup', os.path.basename(vg_path), '-',
                        '-t', str(options.calling_cores)] + options.pileup_opts)
        options.drunner.call(command, work_dir=work_dir, outfile=pu_path_stream)

    # do the calling.
    merged_call_opts = merge_call_opts(path_name, chunk_offset, path_size,
                                       options.call_opts, options.sample_name)
    with open(vcf_path + ".us", "w") as vgcall_stdout, open(vcf_path + ".call_log", "w") as vgcall_stderr:
        command = [['vg', 'call', os.path.basename(vg_path), os.path.basename(pu_path), '-t',
                 str(options.calling_cores)] + str(merged_call_opts).split()]
        options.drunner.call(command, work_dir=work_dir,
                             outfile=vgcall_stdout, errfile=vgcall_stderr)
 
                
def run_vg_genotype(options, xg_path, vg_path, gam_path, path_name, chunk_offset, path_size, work_dir, vcf_path):
                    
    """ Create a VCF with vg genotype """

    # Filter the gam
    # Todo: can we stream directly to indexer?
    filter_gam_path = os.path.join(work_dir, os.path.basename(gam_path) + ".filter")
    with open(filter_gam_path, 'w') as filter_stdout:
        command = ['vg', 'filter', os.path.basename(gam_path), '-t', str(options.calling_cores),
                   '-x', os.path.basename(xg_path)] + options.filter_opts
        options.drunner.call(command, work_dir=work_dir, outfile=filter_stdout)
               
    # Make a gam index
    gam_index_path = os.path.join(work_dir, os.path.basename(filter_gam_path) + ".index")
    command = ['vg', 'index', '-N', os.path.basename(filter_gam_path), '-d', os.path.basename(gam_index_path)]
    options.drunner.call(command, work_dir=work_dir)

    # Do the genotyping
    merged_genotype_opts = merge_call_opts(path_name, chunk_offset, path_size,
                                           options.genotype_options, sample_name, sample_flag = '-s')
    with open(vcf_path + ".us", "w") as vgcall_stdout, open(vcf_path + ".call_log", "w") as vgcall_stderr:

        command=[['vg', 'genotype', os.path.basename(vg_path), os.path.basename(gam_index_path),
                  '-t', str(options.calling_cores), '-v'] + str(merged_genotype_opts).split()]
        options.drunner.call(command, work_dir=work_dir,
                             outfile=vgcall_stdout, errfile=vgcall_stderr)

def call_chunk(job, options, path_name, chunk_i, num_chunks, chunk_offset, clipped_chunk_offset,
               vg_chunk_file_id, gam_chunk_file_id, path_size):
    """ create VCF from a given chunk """
   
    RealTimeLogger.get().info("Running call_chunk on path {} and chunk {}".format(path_name, chunk_i))
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Read our gam files from the store
    vg_path = os.path.join(work_dir, 'chunk_{}_{}.vg'.format(path_name, chunk_offset))
    gam_path = os.path.join(work_dir, 'chunk_{}_{}.gam'.format(path_name, chunk_offset))
    read_from_store(job, options, vg_chunk_file_id, vg_path)
    read_from_store(job, options, gam_chunk_file_id, gam_path)

    # output vcf path
    vcf_path = os.path.join(work_dir, 'chunk_{}_{}.vcf'.format(path_name, chunk_offset))

    # get the xg index, required for vg filter -D.
    # todo? can we avoid this somehow?  is it better to copy in a big xg than regenerate?
    xg_path = os.path.join(work_dir, 'chunk_{}_{}.xg'.format(path_name, chunk_offset))
    options.drunner.call(['vg', 'index', os.path.basename(vg_path), '-x',
                          os.path.basename(xg_path), '-t', str(options.calling_cores)],
                         work_dir = work_dir)
    # Run vg call
    if options.genotype:
        run_vg_genotype(options, xg_path, vg_path, gam_path, path_name, chunk_offset,
                        path_size, work_dir, vcf_path)
    else:
        run_vg_call(options, xg_path, vg_path, gam_path, path_name, chunk_offset,
                    path_size, work_dir, vcf_path)

    # Sort the output
    sort_vcf(options.drunner, vcf_path + ".us", vcf_path)
    options.drunner.call(['rm', vcf_path + '.us'])
    command=['bgzip', '{}'.format(os.path.basename(vcf_path))]
    options.drunner.call(command, work_dir=work_dir)
    command=['tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
    options.drunner.call(command, work_dir=work_dir)
    
    # do the vcf clip
    left_clip = 0 if chunk_i == 0 else options.overlap / 2
    right_clip = 0 if chunk_i == num_chunks - 1 else options.overlap / 2
    clip_path = os.path.join(work_dir, 'chunk_{}_{}_clip.vcf'.format(path_name, chunk_offset))
    with open(clip_path, "w") as clip_path_stream:
        # passing in offset this way pretty hacky, should have its own option
        call_toks = options.call_opts
        offset = 0
        if "-o" in call_toks:
            offset = int(call_toks[call_toks.index("-o") + 1])
        elif "--offset" in call_toks:
            offset = int(call_toks[call_toks.index("--offset") + 1])
        command=['bcftools', 'view', '-r', '{}:{}-{}'.format(
            path_name, offset + clipped_chunk_offset + left_clip + 1,
            offset + clipped_chunk_offset + options.call_chunk_size - right_clip),
                 os.path.basename(vcf_path) + ".gz"]
        options.drunner.call(command, work_dir=work_dir, outfile=clip_path_stream)

    # save clip.vcf files to job store
    clip_file_id = write_to_store(job, options, clip_path)
    
    return clip_file_id

def run_calling(job, options, xg_file_id, alignment_file_id, path_name):
    
    RealTimeLogger.get().info("Running variant calling on path {} from alignment file {}".format(path_name, str(alignment_file_id)))
        
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()

    # Download the input from the store
    xg_path = os.path.join(work_dir, 'graph.vg.xg')
    read_from_store(job, options, xg_file_id, xg_path)
    gam_path = os.path.join(work_dir, '{}.gam'.format(options.sample_name))
    read_from_store(job, options, alignment_file_id, gam_path)

    # Index the gam
    gam_index_path = gam_path + '.index'
    index_cmd = ['vg', 'index', '-N', os.path.basename(gam_path), '-d',
                 os.path.basename(gam_index_path), '-t', str(options.gam_index_cores)]
    options.drunner.call(index_cmd, work_dir=work_dir)

    # Chunk the graph and gam, using the xg and rocksdb indexes
    output_bed_chunks_path = os.path.join(work_dir, 'output_bed_chunks_{}.bed'.format(path_name))
    chunk_cmd = ['vg', 'chunk', '-x', os.path.basename(xg_path),
                 '-a', os.path.basename(gam_index_path), '-c', str(options.chunk_context),
                 '-p', path_name,
                 '-g',
                 '-s', str(options.call_chunk_size),
                 '-o', str(options.overlap),
                 '-b', 'call_chunk_{}'.format(path_name),
                 '-t', str(options.calling_cores),
                 '-R', os.path.basename(output_bed_chunks_path)]
    options.drunner.call(chunk_cmd, work_dir=work_dir)

    # Scrape the BED into memory
    bed_lines = []
    with open(output_bed_chunks_path) as output_bed:
        for line in output_bed:
            toks = line.split('\t')
            if len(toks) > 3:
                bed_lines.append(toks)

    # Infer the size of the path from our BED
    path_size = int(bed_lines[-1][2]) - int(bed_lines[0][1])

    # Go through the BED output of vg chunk, adding a child calling job for
    # each chunk                
    clip_file_ids = []
    for chunk_i, toks in enumerate(bed_lines):
        assert toks[0] == path_name
        chunk_bed_start = int(toks[1])
        chunk_bed_end = int(toks[2])
        chunk_bed_size = chunk_bed_end - chunk_bed_start
        gam_chunk_path = os.path.join(work_dir, os.path.basename(toks[3].strip()))
        assert os.path.splitext(gam_chunk_path)[1] == '.gam'
        vg_chunk_path = os.path.splitext(gam_chunk_path)[0] + '.vg'
        gam_chunk_file_id = write_to_store(job, options, gam_chunk_path)
        vg_chunk_file_id = write_to_store(job, options, vg_chunk_path)
        clipped_chunk_offset = chunk_i * options.call_chunk_size - chunk_i * options.overlap
        
        clip_file_id = job.addChildJobFn(call_chunk, options, path_name, chunk_i, len(bed_lines),
                                         chunk_bed_start, clipped_chunk_offset,
                                         vg_chunk_file_id, gam_chunk_file_id, path_size,
                                         cores=options.calling_cores,
                                         memory=options.calling_mem, disk=options.calling_disk).rv()
        clip_file_ids.append(clip_file_id)


    vcf_gz_tbi_file_id_pair = job.addFollowOnJobFn(merge_vcf_chunks, options, path_name,
                                                   clip_file_ids,
                                                   cores=options.call_chunk_cores,
                                                   memory=options.call_chunk_mem,
                                                   disk=options.call_chunk_disk).rv()
 
    RealTimeLogger.get().info("Completed variant calling on path {} from alignment file {}".format(path_name, str(alignment_file_id)))

    return vcf_gz_tbi_file_id_pair


def merge_vcf_chunks(job, options, path_name, clip_file_ids):
    """ merge a bunch of clipped vcfs created above, taking care to 
    fix up the headers.  everything expected to be sorted already """
    
    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vcf_path = os.path.join(work_dir, path_name + ".vcf")
    
    for chunk_i, clip_file_id in enumerate(clip_file_ids):
        
        # Download clip.vcf file from the store
        clip_path = os.path.join(work_dir, 'clip_{}.vcf'.format(chunk_i))
        read_from_store(job, options, clip_file_id, clip_path)

        if chunk_i == 0:
            # copy everything including the header
            with open(vcf_path, "w") as outfile:
                options.drunner.call(['cat', os.path.basename(clip_path)], outfile=outfile,
                                     work_dir=work_dir)
        else:
            # add on everythin but header
            with open(vcf_path, "a") as outfile:
                options.drunner.call(['bcftools', 'view', '-H', os.path.basename(clip_path)],
                                     outfile=outfile, work_dir=work_dir)

    # add a compressed indexed version
    vcf_gz_file = vcf_path + ".gz"
    with open(vcf_gz_file, "w") as vcf_gz_file_stream:
        command=['bgzip', '-c', '{}'.format(os.path.basename(vcf_path))]
        options.drunner.call(command, work_dir=work_dir, outfile=vcf_gz_file_stream)
    command=['bcftools', 'tabix', '-f', '-p', 'vcf', '{}'.format(os.path.basename(vcf_path+".gz"))]
    options.drunner.call(command, work_dir=work_dir)

    # Save merged vcf files to the job store
    use_out_store = True if options.tool == 'call' else None
    vcf_gz_file_id = write_to_store(job, options, vcf_path+".gz", use_out_store)
    vcf_tbi_file_id = write_to_store(job, options, vcf_path+".gz.tbi", use_out_store)
        
    return vcf_gz_file_id, vcf_tbi_file_id
    
def call_main(options):
    """ no harm in preserving command line access to chunked_call for debugging """
    
    RealTimeLogger.start_master()

    # currently only support (exactly) one path
    require(len(options.chroms) == 1, "Must pass exactly one reference path with --chroms")

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))

    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'call'

    # Throw error if something wrong with IOStore string
    IOStore.get(options.out_store)
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()
    
    with Toil(options) as toil:
        if not toil.options.restart:

            # Upload local files to the remote IO Store
            inputXGFileID = import_to_store(toil, options, options.xg_path)
            inputGamFileID = import_to_store(toil, options, options.gam_path)

            # Make a root job
            root_job = Job.wrapJobFn(run_calling, options, inputXGFileID, inputGamFileID,
                                     options.chroms[0],
                                     cores=options.call_chunk_cores,
                                     memory=options.call_chunk_mem,
                                     disk=options.call_chunk_disk)

            
            # Run the job and store the returned list of output files to download
            toil.start(root_job)
        else:
            toil.restart()
                
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    RealTimeLogger.stop_master()
    
    
