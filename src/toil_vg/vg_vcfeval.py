#!/usr/bin/env python2.7
"""
Thin wrapper for vcfeval.  Expects compressed indexed vcfs as input
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
    parser.add_argument("call_vcf", type=str,
                        help="input vcf (must be bgzipped and have .tbi")
    parser.add_argument("truth_vcf", type=str,
                        help="truth vcf (must be bgzipped and have .tbi")
    parser.add_argument("fasta", type=str,
                        help="fasta file containing the DNA sequence.  Can be"
                        " gzipped with .gz extension")
    parser.add_argument("out_store", type=str,
                        help="vcfeval output store")

    # Add docker options
    add_vcfeval_docker_tool_parse_args(parser)

    # Add common calling options shared with vg_evaluation_pipeline
    vcfeval_parse_args(parser)

    # Add common docker options shared with vg_evaluation pipeline
    add_docker_tool_parse_args(parser)
                        
    return parser.parse_args()

def vcfeval_parse_args(parser):
    """ centralize calling parameters here """

    parser.add_argument("--bed_regions", type=int,
                        help="BED file of regions to consider")
    parser.add_argument("--vcfeval_filtering_opts", type=str,
                        help="Additional filtering options for vcfeval",
                        default="")
    parser.add_argument("--vcfeval_cores", type=int,
                        default=1,
                        help="Cores to use for vcfeval")

def add_vcfeval_docker_tool_parse_args(parser):
    """ vcfeval specific docker options and their defaults """
    parser.add_argument("--rtg_docker", type=list, default=['realtimegenomics/rtg-tools:3.7.1', True],
                        help="dockerfile for rtg_tools")

def get_vcfeval_docker_tool_map(options):
    """ convenience function to parse the above _docker options into a dictionary """

    dmap = dict()
    if not options.no_docker:
        dmap["rtg"] = options.rtg_docker

    return dmap
    
    
def vcfeval(work_dir, call_vcf_name, truth_vcf_name,
            sdf_name, outdir_name, bed_name, options):

    """ create and run the vcfeval command """

    cmd = ['rtg', 'vcfeval', '--calls', call_vcf_name,
           '--baseline', truth_vcf_name,
           '--template', sdf_name, '--output', outdir_name,
           '--threads', str(options.vcfeval_cores)]

    if bed_name is not None:
        cmd += ['--evaluation-regions', os.path.basename(bed_name)]

    if len(options.vcfeval_filtering_opts) > 0:
        cmd += options.vcfeval_filtering_opts.split()

    options.drunner.call(cmd, work_dir=work_dir)

    # get the F1 out of summary.txt
    # expect header on 1st line and data on 3rd
    # todo: be more robust
    with open(os.path.join(work_dir, os.path.basename(outdir_name), "summary.txt")) as sum_file:
        header = sum_file.readline().split()
        line = sum_file.readline()
        data = sum_file.readline().split()
        assert header[-1] == 'F-measure'
        assert len(data) == len(header)
        return float(data[-1])    
    
def run_vcfeval(job, options,
                truth_vcf_name, truth_vcf_id, truth_tbi_id,
                call_vcf_name, call_vcf_id, call_tbi_id,
                fasta_name, fasta_id,
                bed_name, bed_id):                
    """ run vcf_eval, return f1 score """

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    # download into local temp directory
    for id,name in [(truth_vcf_id, truth_vcf_name),
                    (truth_tbi_id, truth_vcf_name + '.tbi'),
                    (call_vcf_id, call_vcf_name),
                    (call_tbi_id, call_vcf_name + '.tbi'),
                    (fasta_id, fasta_name),
                    (bed_id, bed_name)]:
        if id is not None:
            job.fileStore.readGlobalFile(id, os.path.join(work_dir, name))

    # output directory
    out_name = "vcfeval_output"
    # indexed sequence
    sdf_name = fasta_name + ".sdf"
    
    # temporary (i hope) hack:  for some reason rtg wants /data/ on the
    # file paths, while most other images don't seem to care.  add them
    # here for now but would really like to understand why
    if options.drunner.has_tool("rtg"):
        fasta_name = os.path.join('/data', fasta_name)
        sdf_name = os.path.join('/data', sdf_name)
        call_vcf_name = os.path.join('/data', call_vcf_name)
        truth_vcf_name = os.path.join('/data', truth_vcf_name)
        out_name = os.path.join('/data', out_name)
        if bed_name is not None:
            bed_name = os.path.join('/data', bed_name)

    # make an indexed sequence (todo: allow user to pass one in)
    options.drunner.call(['rtg', 'format',  fasta_name, '-o', sdf_name], work_dir=work_dir)    

    # run the vcf_eval command
    f1 = vcfeval(work_dir, call_vcf_name, truth_vcf_name,
                 sdf_name, out_name, bed_name, options)

    # todo : write everything to output store.
    # (just sticking the f1 score in a text file for now)
    f1_path = os.path.join(work_dir, "f1.txt")
    with open(f1_path, "w") as f:
        f.write(str(f1))    
    out_store = IOStore.get(options.out_store)   
    out_store.write_output_file(f1_path, "f1.txt")

    return f1

def main():
    """ command line access to toil vcf eval logic"""
    
    RealTimeLogger.start_master()
    options = parse_args() # This holds the nicely-parsed options object

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = dict(get_docker_tool_map(options).items() +
                               get_vcfeval_docker_tool_map(options).items()))
        
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    # we can relax this down the road by optionally doing some compression/indexing
    assert options.truth_vcf.endswith(".vcf.gz")
    assert options.call_vcf.endswith(".vcf.gz")
    
    with Toil(options) as toil:
        if not toil.options.restart:

            truth_vcf_name = "truth_" + os.path.basename(options.truth_vcf)
            call_vcf_name = "call_" + os.path.basename(options.call_vcf)
            fasta_name = "fa_" + os.path.basename(options.fasta)
            bed_name = "bed_" + os.path.basename(options.bed_regions) if options.bed_regions is not None else None

            # Upload local files to the remote IO Store
            truth_vcf_id = toil.importFile('file://'+options.truth_vcf)
            call_vcf_id = toil.importFile('file://'+options.call_vcf)
            truth_tbi_id = toil.importFile('file://'+options.truth_vcf + '.tbi')
            call_tbi_id = toil.importFile('file://'+options.call_vcf + '.tbi')            
            fasta_id = toil.importFile('file://'+options.fasta)
            bed_id = toil.importFile('file://'+options.bed_regions) if options.bed_regions is not None else None
            
            # Make a root job
            root_job = Job.wrapJobFn(run_vcfeval, options,
                                     truth_vcf_name, truth_vcf_id, truth_tbi_id,
                                     call_vcf_name, call_vcf_id, call_tbi_id,
                                     fasta_name, fasta_id, bed_name, bed_id,
                                     cores=2, memory="5G", disk="2G")
            
            # Run the job
            f1 = toil.start(root_job)
        else:
            f1 = toil.restart()

        print("F1 Score : {}".format(f1))
                
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
