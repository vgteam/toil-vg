#!/usr/bin/env python2.7
"""
Thin wrapper for vcfeval, as a convenience to stick a vcfeval output directory
along with the other toil-vg output.  Can be run standalone as well.
"""
from __future__ import print_function
import argparse, sys, os, os.path, random, subprocess, shutil, itertools, glob
import json, time, timeit, errno
from uuid import uuid4
import logging

from toil.common import Toil
from toil.job import Job
from toil_vg.vg_common import *

logger = logging.getLogger(__name__)

def vcfeval_subparser(parser):
    """
    Create a subparser for vcf evaluation.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)

    # General options
    parser.add_argument("call_vcf", type=str,
                        help="input vcf (must be bgzipped and have .tbi")
    parser.add_argument("vcfeval_baseline", type=str,
                        help="truth vcf (must be bgzipped and have .tbi")
    parser.add_argument("vcfeval_fasta", type=str,
                        help="fasta file containing the DNA sequence.  Can be"
                        " gzipped with .gz extension")
    parser.add_argument("out_store",
                            help="output store.  All output written here. Path specified using same syntax as toil jobStore")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common calling options shared with toil_vg pipeline
    vcfeval_parse_args(parser)

    # Add common docker options shared with toil_vg pipeline
    add_container_tool_parse_args(parser)

def vcfeval_parse_args(parser):
    """ centralize calling parameters here """

    parser.add_argument("--vcfeval_bed_regions", type=str,
                        help="BED file of regions to consider")
    parser.add_argument("--vcfeval_opts", type=str,
                        help="Additional options for vcfeval (wrapped in \"\")",
                        default=None)
    parser.add_argument("--vcfeval_cores", type=int,
                        default=1,
                        help="Cores to use for vcfeval")
    
def vcfeval(job, work_dir, call_vcf_name, vcfeval_baseline_name,
            sdf_name, outdir_name, bed_name, options):

    """ create and run the vcfeval command """

    cmd = ['rtg', 'vcfeval', '--calls', call_vcf_name,
           '--baseline', vcfeval_baseline_name,
           '--template', sdf_name, '--output', outdir_name,
           '--threads', str(options.vcfeval_cores),
           '--vcf-score-field', 'QUAL']

    if bed_name is not None:
        cmd += ['--evaluation-regions', bed_name]

    if options.vcfeval_opts:
        cmd += options.vcfeval_opts

    options.drunner.call(job, cmd, work_dir=work_dir)

    # get the F1 out of summary.txt
    # expect header on 1st line and data on 3rd and below
    # we take the best F1 found over these lines (which should correspond to best
    # point on quality ROC curve)
    # todo: be more robust
    f1 = None
    with open(os.path.join(work_dir, os.path.basename(outdir_name), "summary.txt")) as sum_file:
        header = sum_file.readline().split()
        assert header[-1] == 'F-measure'        
        line = sum_file.readline()
        for line in sum_file:
            data = line.strip().split()
            assert len(data) == len(header)
            line_f1 = float(data[-1])
            if f1 is None or line_f1 > f1:
                f1 = line_f1
    return f1
    
def run_vcfeval(job, options, vcf_tbi_id_pair, vcfeval_baseline_id, vcfeval_baseline_tbi_id, 
                fasta_id, bed_id):                
    """ run vcf_eval, return f1 score """

    # make a local work directory
    work_dir = job.fileStore.getLocalTempDir()

    call_vcf_id, call_tbi_id = vcf_tbi_id_pair[0], vcf_tbi_id_pair[1]

    vcfeval_baseline_name = "truth_" + os.path.basename(options.vcfeval_baseline)
    call_vcf_name = "calls.vcf.gz"
    fasta_name = "fa_" + os.path.basename(options.vcfeval_fasta)
    bed_name = "bed_" + os.path.basename(options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None   
    # download into local temp directory
    for id,name in [(vcfeval_baseline_id, vcfeval_baseline_name),
                    (vcfeval_baseline_tbi_id, vcfeval_baseline_name + '.tbi'),
                    (call_vcf_id, call_vcf_name),
                    (call_tbi_id, call_vcf_name + '.tbi'),
                    (fasta_id, fasta_name),
                    (bed_id, bed_name)]:
        if id is not None:
            read_from_store(job, options, id, os.path.join(work_dir, name))

    try:
        out_tag = '{}_vcfeval_output'.format(options.sample_name)
    except:
        out_tag = 'vcfeval_output'
        
    # output directory
    out_name = out_tag
    # indexed sequence
    sdf_name = fasta_name + ".sdf"
    
    # make an indexed sequence (todo: allow user to pass one in)
    options.drunner.call(job, ['rtg', 'format',  fasta_name, '-o', sdf_name], work_dir=work_dir)    

    # run the vcf_eval command
    f1 = vcfeval(job, work_dir, call_vcf_name, vcfeval_baseline_name,
                 sdf_name, out_name, bed_name, options)

    # copy results to the output store
    # 1) vcfeval_output_f1.txt (used currently by tests script)
    f1_path = os.path.join(work_dir, "f1.txt")    
    with open(f1_path, "w") as f:
        f.write(str(f1))
    write_to_store(job, options, os.path.join(work_dir, out_tag, 'summary.txt'), use_out_store = True,
                   out_store_key = '{}_summary.txt'.format(out_tag))
    # 2) vcfeval_output_summary.txt
    write_to_store(job, options, f1_path, use_out_store = True,
                   out_store_key = '{}_f1.txt'.format(out_tag))
    options.drunner.call(job, ['tar', 'czf', out_tag + '.tar.gz', out_tag], work_dir = work_dir)
    # 3) vcfeval_output.tar.gz -- whole shebang
    write_to_store(job, options, os.path.join(work_dir, out_tag + '.tar.gz'), use_out_store = True)
    # 4) truth VCF
    write_to_store(job, options, os.path.join(work_dir, vcfeval_baseline_name), use_out_store = True)

    return f1

def vcfeval_main(options):
    """ command line access to toil vcf eval logic"""
    
    # make the docker runner
    options.drunner = ContainerRunner(
        container_tool_map = get_container_tool_map(options))
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    # we can relax this down the road by optionally doing some compression/indexing
    assert options.vcfeval_baseline.endswith(".vcf.gz")
    assert options.call_vcf.endswith(".vcf.gz")
    
    with Toil(options) as toil:
        if not toil.options.restart:

            start_time = timeit.default_timer()
            
            # Upload local files to the remote IO Store
            vcfeval_baseline_id = import_to_store(toil, options, options.vcfeval_baseline)
            call_vcf_id = import_to_store(toil, options, options.call_vcf)
            vcfeval_baseline_tbi_id = import_to_store(toil, options, options.vcfeval_baseline + '.tbi')
            call_tbi_id = import_to_store(toil, options, options.call_vcf + '.tbi')            
            fasta_id = import_to_store(toil, options, options.vcfeval_fasta)
            bed_id = import_to_store(toil, options, options.vcfeval_bed_regions) if options.vcfeval_bed_regions is not None else None

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_vcfeval, options,
                                     (call_vcf_id, call_tbi_id),
                                     vcfeval_baseline_id, vcfeval_baseline_tbi_id,
                                     fasta_id, bed_id,
                                     cores=options.vcfeval_cores, memory=options.vcfeval_mem,
                                     disk=options.vcfeval_disk)

            # Run the job
            f1 = toil.start(root_job)
        else:
            f1 = toil.restart()

        print("F1 Score : {}".format(f1))
                
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
    
