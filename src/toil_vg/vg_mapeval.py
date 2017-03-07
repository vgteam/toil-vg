#!/usr/bin/env python2.7
"""
vg_mapeval.py: Compare alignment positions from gam or bam to a truth set
that was created with vg sim --gam

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

def mapeval_subparser(parser):
    """
    Create a subparser for mapeval.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # General options
    parser.add_argument("out_store",
                        help="output store.  All output written here. Path specified using same syntax as toil jobStore")
    parser.add_argument("xg",
                        help="xg corresponding to gam")
    parser.add_argument("gam",
                        help="aligned reads to compare to truth")
    parser.add_argument("reads_gam", default=None,
                        help="reads in GAM format as output by toil-vg sim")
    parser.add_argument("truth", default=None,
                        help="list of true positions of reads as output by toil-vg sim")

    parser.add_argument("--mapeval_threshold", type=int, default=100,
                        help="distance between alignment and true position to be called correct")

    parser.add_argument("--bwa", action="store_true",
                        help="run bwa mem on the reads, and add to comparison")
    parser.add_argument("--bwa-paired", action="store_true",
                        help="run bwa mem paired end as well")
    parser.add_argument('--fasta', default=None,
                        help="fasta sequence file (required for bwa)")
    parser.add_argument('--gam_reads', default=None,
                        help="reads in GAM format (required for bwa)")

    parser.add_argument("--bwa_opts", type=str,
                        help="arguments for bwa mem (wrapped in \"\").")
        

    # Add common options shared with everybody
    add_common_vg_parse_args(parser)

    # Add common docker options
    add_docker_tool_parse_args(parser)

def run_bwa_index(job, options, gam_file_id, fasta_file_id, bwa_index_ids):
    """
    Make a bwa index for a fast sequence if not given in input. then run bwa mem
    """
    if not bwa_index_ids:
        bwa_index_ids = dict()
        work_dir = job.fileStore.getLocalTempDir()
        fasta_file = os.path.join(work_dir, os.path.basename(options.fasta))
        read_from_store(job, options, fasta_file_id, fasta_file)
        cmd = ['bwa', 'index', os.path.basename(fasta_file)]
        options.drunner.call(job, cmd, work_dir = work_dir)
        for idx_file in glob.glob('{}.*'.format(fasta_file)):
            bwa_index_ids[idx_file[len(fasta_file):]] = write_to_store(job, options, idx_file)

    bwa_pos_file_id = None
    bwa_pair_pos_file_id = None
                    
    if options.bwa:
        bwa_pos_file_id = job.addChildJobFn(run_bwa_mem, options, gam_file_id, bwa_index_ids, False,
                                            cores=options.alignment_cores, memory=options.alignment_mem,
                                            disk=options.alignment_disk).rv()
    if options.bwa_paired:
        bwa_pair_pos_file_id = job.addChildJobFn(run_bwa_mem, options, gam_file_id, bwa_index_ids, True,
                                                 cores=options.alignment_cores, memory=options.alignment_mem,
                                                 disk=options.alignment_disk).rv()

    return bwa_pos_file_id, bwa_pair_pos_file_id

    
def run_bwa_mem(job, options, gam_file_id, bwa_index_ids, paired_mode):
    """ run bwa-mem on reads in a gam.  optionally run in paired mode
    return id of positions file
    """

    work_dir = job.fileStore.getLocalTempDir()

    # read the gam file
    gam_file = os.path.join(work_dir, os.path.basename(options.reads_gam))
    read_from_store(job, options, gam_file_id, gam_file)

    # and the index files
    fasta_file = os.path.join(work_dir, os.path.basename(options.fasta))
    for suf, idx_id in bwa_index_ids.items():
        RealtimeLogger.info("reading index {}".format('{}{}'.format(fasta_file, suf)))
        read_from_store(job, options, idx_id, '{}{}'.format(fasta_file, suf))

    # output positions file
    positions_file = os.path.join(work_dir, 'bwa.pos')
    
    # if we're paired, must make some split files
    if paired_mode:

        # convert to json (todo: have docker image that can do vg and jq)
        json_file = gam_file + '.json'
        cmd = ['vg', 'view', '-a', os.path.basename(gam_file)]
        with open(json_file, 'w') as out_json:
            options.drunner.call(job, cmd, work_dir = work_dir, outfile = out_json)

        sim_fq_files = [None, os.path.join(work_dir, 'sim_1.fq.gz'),
                        os.path.join(work_dir, 'sim_2.fq.gz')]

        # jq docker image is another that requires the /data/.  really need to figure
        # out more general approach
        json_path = os.path.basename(json_file)
        if options.drunner.has_tool("jq"):
            json_path = os.path.join('/data', json_path)
        
        # make a fastq for each end of pair
        for i in [1, 2]:
            # extract paired end with jq
            cmd = ['jq', '-cr', 'select(.name | test("_{}$"))'.format(i), json_path]
            end_file = json_file + '.{}'.format(i)
            with open(end_file, 'w') as end_out:
                options.drunner.call(job, cmd, work_dir = work_dir, outfile = end_out)

            cmd = [['vg', 'view', '-JaG', os.path.basename(end_file)]]
            cmd.append(['vg', 'view', '-X', '-'])
            cmd.append(['sed', 's/_{}$//'.format(i)])
            cmd.append(['gzip'])

            with open(sim_fq_files[i], 'w') as sim_out:
                options.drunner.call(job, cmd, work_dir = work_dir, outfile = sim_out)

            os.remove(end_file)

        # run bwa-mem on the paired end input
        cmd = [['bwa', 'mem', '-t', str(options.alignment_cores), os.path.basename(fasta_file),
                os.path.basename(sim_fq_files[1]), os.path.basename(sim_fq_files[2])] + options.bwa_opts]
        cmd.append(['grep', '-v', '^@'])
        cmd.append(['perl', '-ne', '@val = split("\t", $_); print @val[0] . "_" . (@val[1] & 64 ? "1" : @val[1] & 128 ? "2" : "?"), "\t" . @val[2] . "\t" . @val[3] . "\t" . @val[4] . "\n";'])
        cmd.append(['sort'])
        
        with open(positions_file, 'w') as out_pos:
            options.drunner.call(job, cmd, work_dir = work_dir, outfile = out_pos)

    # single end
    else:

        # extract reads from gam.  as above, need to have single docker container (which shouldn't be
        # a big deal) to run all these chained command and avoid huge files on disk
        extracted_reads_file = os.path.join(work_dir, 'extracted_reads')
        cmd = ['vg', 'view', '-X', os.path.basename(gam_file)]
        with open(extracted_reads_file, 'w') as out_ext:
            options.drunner.call(job, cmd, work_dir = work_dir, outfile = out_ext)

        # run bwa-mem on single end input
        cmd = [['bwa', 'mem', '-t', str(options.alignment_cores), os.path.basename(fasta_file),
                os.path.basename(extracted_reads_file)] + options.bwa_opts]
        cmd.append(['grep', '-v', '^@'])
        cmd.append(['cut', '-f', '1,3,4,5'])
        cmd.append(['sort'])
        with open(positions_file, 'w') as out_pos:
            options.drunner.call(job, cmd, work_dir = work_dir, outfile = out_pos)

    # return our id for the output positions file
    pos_file_id = write_to_store(job, options, positions_file)
    return pos_file_id

def get_gam_positions(job, options, work_dir, xg_file, gam_file, out_pos_file):
    """
    extract positions from gam, return id of positions file
    (lots of duplicated code with vg_sim, should merge?)

    """

    # go through intermediate json file until docker worked out
    gam_annot_json = gam_file + '.json'
    cmd = [['vg', 'annotate', '-p', '-a', os.path.basename(gam_file), '-x', os.path.basename(xg_file)]]
    cmd.append(['vg', 'view', '-aj', '-'])
    with open(gam_annot_json, 'w') as output_annot_json:
        options.drunner.call(job, cmd, work_dir = work_dir, outfile=output_annot_json)

    # jq docker image is another that requires the /data/.  really need to figure
    # out more general approach
    json_path = os.path.basename(gam_annot_json)
    if options.drunner.has_tool("jq"):
        json_path = os.path.join('/data', json_path)
        
    # turn the annotated gam json into truth positions, as separate command since
    # we're going to use a different docker container.  (Note, would be nice to
    # avoid writing the json to disk)        
    jq_cmd = [['jq', '-c', '-r', '[.name, .refpos[0].name, .refpos[0].offset,'
               'if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv',
               gam_annot_json]]
    jq_cmd.append(['sed', 's/null/0/g'])

    with open(out_pos_file + '.unsorted', 'w') as out_pos:
        options.drunner.call(job, jq_cmd, work_dir = work_dir, outfile=out_pos)

    # get rid of that big json asap
    os.remove(gam_annot_json)

    # sort the positions file (not piping due to memory fears)
    sort_cmd = ['sort', os.path.basename(out_pos_file) + '.unsorted']
    with open(out_pos_file, 'w') as out_pos:
        options.drunner.call(job, sort_cmd, work_dir = work_dir, outfile = out_pos)
    

def compare_positions(job, options, truth_file_id, pos_file_id):
    """
    this is essentially pos_compare.py from vg/scripts
    return output file id.
    """
    work_dir = job.fileStore.getLocalTempDir()

    true_pos_file = os.path.join(work_dir, 'true.pos')
    read_from_store(job, options, truth_file_id, true_pos_file)
    test_pos_file = os.path.join(work_dir, 'test_pos')
    read_from_store(job, options, pos_file_id, test_pos_file)

    out_file = os.path.join(work_dir, 'out.compare')

    join_file = os.path.join(work_dir, 'join_pos')
    cmd = ['join', os.path.basename(true_pos_file), os.path.basename(test_pos_file)]
    

    with open(true_pos_file) as truth, open(test_pos_file) as test, \
         open(out_file, 'w') as out:
        for truth_line, test_line in zip(truth, test):
            true_fields = truth_line.split()
            test_fields = test_line.split()
            # every input has a true position, and if it has less than the expected number of fields we assume alignment failed
            true_read_name = true_fields[0]
            if len(true_fields) + len(test_fields) != 7:
                out.write('{}, 0, 0\n'.format(true_read_name))
                continue
            
            true_chr = true_fields[1]
            true_pos = int(true_fields[2])
            aln_read_name = test_fields[0]
            assert aln_read_name == true_read_name
            aln_chr = test_fields[1]
            aln_pos = int(test_fields[2])
            aln_mapq = int(test_fields[3])
            aln_correct = 1 if aln_chr == true_chr and abs(true_pos - aln_pos) < options.mapeval_threshold else 0

            out.write('{}, {}, {}\n'.format(aln_read_name, aln_correct, aln_mapq))
        
        # make sure same length
        try:
            iter(true_pos_file).next()
            raise RuntimeError('position files have different lengths')
        except:
            pass
        try:
            iter(test_pos_file).next()
            raise RuntimeError('position files have different lengths')
        except:
            pass
        
    out_file_id = write_to_store(job, options, out_file)
    return out_file_id

def run_map_eval(job, options, xg_file_id, gam_file_id, reads_gam_file_id,
                 fasta_file_id, bwa_file_ids, true_pos_file_id):
    """ run the mapping comparison.  Dump some tables into the outstore """
    
    # run bwa if requested
    bwa_pos_ids = None
    if options.bwa or options.bwa_paired:
        bwa_pos_ids = job.addChildJobFn(run_bwa_index, options, reads_gam_file_id,
                                        fasta_file_id, bwa_file_ids,
                                        cores=options.alignment_cores, memory=options.alignment_mem,
                                        disk=options.alignment_disk).rv()

    # get the gam positions
    work_dir = job.fileStore.getLocalTempDir()

    # read the gam file
    gam_file = os.path.join(work_dir, os.path.basename(options.gam))
    read_from_store(job, options, gam_file_id, gam_file)

    # and the xg
    xg_file = os.path.join(work_dir, os.path.basename(options.xg))
    read_from_store(job, options, xg_file_id, xg_file)

    # extract positions
    gam_pos_file = gam_file + '.pos'
    get_gam_positions(job, options, work_dir, xg_file, gam_file, gam_pos_file)
    gam_pos_file_id = write_to_store(job, options, gam_pos_file)

    # compare all our positions
    comparison_results = job.addFollowOnJobFn(run_map_eval_compare, options, true_pos_file_id,
                                              gam_pos_file_id, bwa_pos_ids,
                                              cores=options.misc_cores, memory=options.misc_mem,
                                              disk=options.misc_disk).rv()

    return comparison_results

def run_map_eval_compare(job, options, true_pos_file_id, gam_pos_file_id,
                         bwa_pos_ids):
    """
    run compare on the positions
    """

    gam_compare_id = job.addChildJobFn(compare_positions, options, true_pos_file_id, gam_pos_file_id,
                                       cores=options.misc_cores, memory=options.misc_mem,
                                       disk=options.misc_disk).rv()

    if bwa_pos_ids:
        bwa_pos_file_id, bwa_pair_pos_file_id = bwa_pos_ids[0], bwa_pos_ids[1]
    else:
        bwa_pos_file_id, bwa_pair_pos_file_id = None, None
    
    bwa_compare_id = None
    if bwa_pos_file_id:
        bwa_compare_id = job.addChildJobFn(compare_positions, options, true_pos_file_id, bwa_pos_file_id,
                                       cores=options.misc_cores, memory=options.misc_mem,
                                       disk=options.misc_disk).rv()

    bwa_pair_compare_id = None
    if bwa_pair_pos_file_id:
        bwa_pair_compare_id = job.addChildJobFn(compare_positions, options, true_pos_file_id, bwa_pair_pos_file_id,
                                                cores=options.misc_cores, memory=options.misc_mem,
                                                disk=options.misc_disk).rv()

    return job.addFollowOnJobFn(run_process_comparisons, options, gam_compare_id,
                                bwa_compare_id, bwa_pair_compare_id,
                                cores=options.misc_cores, memory=options.misc_mem,
                                disk=options.misc_disk).rv()

def run_process_comparisons(job, options, gam_compare_id, bwa_compare_id, bwa_pair_compare_id):
    """
    do something with all the comparison tables.  (just writing them to the outstore for now)
    """

    work_dir = job.fileStore.getLocalTempDir()

    # make the results.tsv
    results_file = os.path.join(work_dir, 'results.tsv')
    with open(results_file, 'w') as out_results:
        out_results.write('correct\tmq\taligner\n')

        def write_tsv(comp_file, a):
            with open(comp_file) as comp_in:
                for line in comp_in:
                    toks = line.rstrip().split(', ')
                    out_results.write('{}\t{}\t"{}"\n'.format(toks[1], toks[2], a))

        # dump our compare files out
        gam_compare_file = os.path.join(work_dir, os.path.basename(options.gam) + '.compare')
        read_from_store(job, options, gam_compare_id, gam_compare_file)
        write_to_store(job, options, gam_compare_file, use_out_store = True)
        write_tsv(gam_compare_file, 'vg')
        
        if bwa_compare_id:
            bwa_compare_file = os.path.join(work_dir, 'bwa.compare')
            read_from_store(job, options, bwa_compare_id, bwa_compare_file)
            write_to_store(job, options, bwa_compare_file, use_out_store = True)
            write_tsv(bwa_compare_file, 'bwa')

        if bwa_pair_compare_id:
            bwa_pair_compare_file = os.path.join(work_dir, 'bwa_pe.compare')
            read_from_store(job, options, bwa_pair_compare_id, bwa_pair_compare_file)
            write_to_store(job, options, bwa_pair_compare_file, use_out_store = True)
            write_tsv(bwa_pair_compare_file, 'bwa-pe')

    write_to_store(job, options, results_file, use_out_store = True)

def mapeval_main(options):
    """
    Wrapper for vg map. 
    """

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))

    if options.bwa or options.bwa_paired:
        require(options.reads_gam, '--reads_gam required for bwa')
        require(options.fasta, '--fasta required for bwa')

    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'mapeval'

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
            inputXGFileID = import_to_store(toil, options, options.xg)
            inputGAMFileID = import_to_store(toil, options, options.gam)
            if options.reads_gam:
                inputReadsGAMFileID = import_to_store(toil, options, options.reads_gam)
            else:
                inputReadsGAMFileID = None
            if options.fasta:
                inputBwaIndexIDs = dict()
                for suf in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
                    fidx = '{}{}'.format(options.fasta, suf)
                    if os.path.exists(fidx):
                        inputBwaIndexIDs[suf] = import_to_store(toil, options, fidx)
                    else:
                        inputBwaIndexIDs = None
                        break
                if not inputBwaIndexIDs:
                    inputFastaID = import_to_store(toil, options, options.fasta)
            else:
                inputFastaID = None
                inputBwaIndexIDs = None
            inputTruePosFileID = import_to_store(toil, options, options.truth)

            end_time = timeit.default_timer()
            logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))

            # Make a root job
            root_job = Job.wrapJobFn(run_map_eval, options, inputXGFileID,
                                     inputGAMFileID,
                                     inputReadsGAMFileID,
                                     inputFastaID,
                                     inputBwaIndexIDs,
                                     inputTruePosFileID,
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
    
