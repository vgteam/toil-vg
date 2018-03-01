#!/usr/bin/env python2.7
"""
vg_mapeval.py: Compare alignment positions from gam or bam to a truth set
that was created with vg sim --gam

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string, math
import urlparse
import getpass
import pdb
import gzip
import logging
import copy
from collections import Counter

from math import ceil
from subprocess import Popen, PIPE

try:
    import numpy as np
    from sklearn.metrics import roc_auc_score, average_precision_score, r2_score, roc_curve
    have_sklearn = True
except:
    have_sklearn = False

import tsv

from toil.common import Toil
from toil.job import Job
from toil.realtimeLogger import RealtimeLogger
from toil_vg.vg_common import require, make_url, remove_ext,\
    add_common_vg_parse_args, add_container_tool_parse_args, get_vg_script
from toil_vg.vg_map import map_parse_args, run_mapping
from toil_vg.vg_index import run_indexing
from toil_vg.context import Context, run_write_info_to_outstore

logger = logging.getLogger(__name__)

def mapeval_subparser(parser):
    """
    Create a subparser for mapeval.  Should pass in results of subparsers.add_parser()
    """

    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser)
    
    # Add the out_store
    # TODO: do this at a higher level?
    # Or roll into Context?
    parser.add_argument('out_store',
                        help='output store.  All output written here. Path specified using same syntax as toil jobStore')
    
    # Add mapeval stuff
    add_mapeval_options(parser)
    
    # Add common docker options
    add_container_tool_parse_args(parser)
    
def add_mapeval_options(parser):
    """
    Add the mapeval options to the given argparse parser.
    """
    
    # General options
    parser.add_argument('--truth', type=make_url, default=None,
                        help='list of true positions of reads as output by toil-vg sim'
                        ' (by default positions extracted from --gam_input_reads or --bam_input_reads)')        
    parser.add_argument('--gams', nargs='+', type=make_url, default=[],
                        help='aligned reads to compare to truth.  specify xg index locations with --index-bases')
    parser.add_argument("--index-bases", nargs='+', type=make_url, default=[],
                        help='use in place of gams to perform alignment.  will expect '
                        '<index-base>.gcsa, <index-base>.lcp and <index-base>.xg to exist')
    parser.add_argument('--use-gbwt', action='store_true',
                        help='also import <index-base>.gbwt and use it during alignment')
    parser.add_argument('--use-snarls', action='store_true',
                        help='also import <index-base>.snarls and use it during multipath alignment')
    parser.add_argument('--vg-graphs', nargs='+', type=make_url, default=[],
                        help='vg graphs to use in place of gams or indexes.  indexes'
                        ' will be built as required')
    parser.add_argument('--gam-names', nargs='+', default=[],
                        help='a name for each gam passed in --gams/graphs/index-bases')
    parser.add_argument('--bams', nargs='+', type=make_url, default=[],
                        help='aligned reads to compare to truth in BAM format')
    parser.add_argument('--bam-names', nargs='+', default=[],
                        help='a name for each bam passed in with --bams')
    parser.add_argument('--pe-bams', nargs='+', type=make_url, default=[],
                        help='paired end aligned reads t compare to truth in BAM format')
    parser.add_argument('--pe-bam-names', nargs='+', default=[],
                        help='a name for each bam passed in with --pe-bams')
    parser.add_argument('--paired-only', action='store_true',
                        help='only do paired-end alignment (default is to do single and paired)')
    parser.add_argument('--single-only', action='store_true',
                        help='only do single-end alignment (default is to do single and paired)')
    
    parser.add_argument('--mapeval-threshold', type=int, default=200,
                        help='distance between alignment and true position to be called correct')

    parser.add_argument('--bwa', action='store_true',
                        help='run bwa mem on the reads, and add to comparison')
    parser.add_argument('--fasta', type=make_url, default=None,
                        help='fasta sequence file (required for bwa. if a bwa index exists for this file, it will be used)')
    parser.add_argument('--bwa-opts', type=str,
                        help='arguments for bwa mem (wrapped in \"\").')
    
    # We can compare all the scores against those from a particular GAM, if asked.
    parser.add_argument('--compare-gam-scores', default=None,
                        help='compare scores against those in the given named GAM')

    parser.add_argument('--ignore-quals', action='store_true',
                        help='never use quality adjusted alignment. ' 
                        'necessary if using --multipath on reads not from trained simulator')

    parser.add_argument('--multipath-only', action='store_true',
                        help='run only mpmap and not map (--multipath will run both and using neither will just run map)')

    parser.add_argument('--more-mpmap-opts', nargs='+', default=[],
                        help='additional batches of mpmap options to try')

    parser.add_argument('--gam-input-xg', type=make_url, default=None,
                        help= 'If extracting truth positions from --input_gam_reads, specify corresponding xg for annotation')
                        
    # We also need to have these options to make lower-level toil-vg code happy
    # with the options namespace we hand it.
    
    # Add mapping options
    map_parse_args(parser)
    
    # Add common options shared with everybody
    add_common_vg_parse_args(parser)
    
def get_default_mapeval_options():
    """
    Return an argparse Namespace populated with the default mapeval option
    values.
    
    Requires the required positional truth file argument.
    
    Can be modified and then passed to make_mapeval_plan(), so you can use
    mapeval as part of a larger program.
    
    """
    
    # Make a parser
    parser = argparse.ArgumentParser()
    # Stick our arguments on it
    add_mapeval_options(parser)
    # And parse nothing but mandatory arguments
    return parser.parse_args([])
    
def validate_options(options):
    """
    Throw an error if an invalid combination of options has been selected.
    """

    # need to have input reads coming from somewhere
    require(sum(map(lambda x : x is not None,
                    [options.gam_input_reads, options.bam_input_reads, options.fastq])) == 1,
            'one of --gam_input_reads or --fastq or --bam_input_reads required for input')

    # annotation is not an option when reading fastq
    require(not options.fastq or options.truth,
            '--truth required with --fastq input')

    # only one or two fastqs accepted
    require(not options.fastq or len(options.fastq) in [1,2],
            'only 1 or two fastqs accepted with --fatsq')

    # only gzipped fastqs accpeted
    require(not options.fastq or all(map(lambda x : x.endswith('.gz'), options.fastq)),
            'only gzipped fastqs (ending with .gz) accepted by --fastq')
            
    # check bwa / bam input parameters.  
    if options.bwa:
        require(options.fasta, '--fasta required for bwa')
    if options.bams:
        require(options.bam_names and len(options.bams) == len(options.bam_names),
                 '--bams and --bam-names must have same number of inputs')
    if options.pe_bams:
        require(options.pe_bam_names and len(options.pe_bams) == len(options.pe_bam_names),
                 '--pe-bams and --pe-bam-names must have same number of inputs')

    # some options from toil-vg map are disabled on the command line
    # this can be eventually cleaned up a bit better 
    require(not options.interleaved,
            '--interleaved disabled in toil-vg mapeval; a single --fastq is always assumed interleaved and two are always assumed paired')

    # accept graphs or indexes in place of gams
    require(options.gams or options.index_bases or options.vg_graphs,
            'one of --vg-graphs, --index-bases or --gams must be used to specifiy vg input')

    if options.use_gbwt:
        require(not options.gams,
                '--use-gbwt cannot be used with pre-aligned GAMs in --gams')
                
    if options.use_snarls:
        require(options.multipath or options.multipath_only,
                '--use-snarls only affects the multipath mapper (--multipath or --multipath-only)')

    if options.gams:
        require(len(options.index_bases) == len(options.gams),
                '--index-bases must be used along with --gams to specify xg locations')
    if options.vg_graphs:
        require(not options.gams and not options.index_bases,
                'if --vg-graphs specified, --gams and --index-bases must not be used')

    # must have a name for each graph/index/gam
    if options.gams:
        require(options.gam_names and len(options.gams) == len(options.gam_names),
                 '--gams and --gam_names must have same number of inputs')
    if options.vg_graphs:
        require(options.gam_names and len(options.vg_graphs) == len(options.gam_names),
                 '--vg-graphs and --gam_names must have same number of inputs')
    if options.index_bases:
        require(options.gam_names and len(options.index_bases) == len(options.gam_names),
                 '--index-bases and --gam_names must have same number of inputs')
                 
    # Make sure names are unique so we can use them as dict keys and file names
    names = []
    if options.bam_names:
        names += options.bam_names
    if options.pe_bam_names:
        names += options.pe_bam_names
    if options.gam_names:
        names += options.gam_names
    require(len(names) == len(set(names)), 'all names must be unique')

    require(options.gam_input_reads is None or options.bam_input_reads is None,
            '--gam_input_reads and --bam_input_reads cannot both be specified')

    require(options.truth or options.bam_input_reads or options.gam_input_xg,
            '--gam-input-xg must be used to specify xg index to annotate --gam_input_reads')
    
def parse_int(value):
    """
    Parse an int, interpreting an empty string as 0.
    """
    
    return int(value) if value.strip() != '' else 0

def run_bwa_index(job, context, fasta_file_id, bwa_index_ids):
    """
    Make a bwa index for a fast sequence if not given in input. then run bwa mem
    
    Retuns a pair of BAM IDs, one for unpaired and one for paired. Either will be None if the corresponding flag is false.
    """
    if not bwa_index_ids:
        bwa_index_ids = dict()
        work_dir = job.fileStore.getLocalTempDir()
        # Download the FASTA file to be indexed
        # It would be nice to name it the same as the actual input FASTA but we'd have to peek at the options
        fasta_file = os.path.join(work_dir, 'toindex.fa')
        job.fileStore.readGlobalFile(fasta_file_id, fasta_file)
        cmd = ['bwa', 'index', os.path.basename(fasta_file)]
        context.runner.call(job, cmd, work_dir = work_dir)
        for idx_file in glob.glob('{}.*'.format(fasta_file)):
            # Upload all the index files created, and store their IDs under their extensions
            bwa_index_ids[idx_file[len(fasta_file):]] = context.write_intermediate_file(job, idx_file)

    return bwa_index_ids

def run_bam_to_fastq(job, context, bam_file_id, paired_mode, add_paired_suffix=False):
    """
    convert a bam to fastq (or pair of fastqs).  add_suffix will stick a _1 or _2 on
    paired reads (needed for vg, but not bwa)
    """
    work_dir = job.fileStore.getLocalTempDir()

    # read the bam file
    bam_file = os.path.join(work_dir, 'input.bam')
    job.fileStore.readGlobalFile(bam_file_id, bam_file)
    
    # if we're paired, must make some split files
    if paired_mode:
        sim_fq_files = [os.path.join(work_dir, 'sim_1{}.fq'.format('s' if add_paired_suffix else '')),
                        os.path.join(work_dir, 'sim_2{}.fq'.format('s' if add_paired_suffix else ''))]
        cmd = ['samtools', 'fastq', os.path.basename(bam_file),
               '-1', os.path.basename(sim_fq_files[0]),
               '-2', os.path.basename(sim_fq_files[1])]
        if add_paired_suffix:
            cmd += ['-N']
        else:
            cmd += ['-n']
        context.runner.call(job, cmd, work_dir = work_dir)
        # we change /1 /2 --> _1 _2 to be compatible with rest of mapeval
        gzip_cmd = [['sed', os.path.basename(sim_fq_files[0]), '-e', 's/\/1/_1/g'], ['gzip', '-c']]
        with open(sim_fq_files[0] + '.gz', 'w') as gz_file:
            context.runner.call(job, gzip_cmd, work_dir = work_dir, outfile = gz_file)
        gzip_cmd = [['sed', os.path.basename(sim_fq_files[1]), '-e', 's/\/2/_2/g'], ['gzip', '-c']]
        with open(sim_fq_files[1] + '.gz', 'w') as gz_file:
            context.runner.call(job, gzip_cmd, work_dir = work_dir, outfile = gz_file)
        return [context.write_intermediate_file(job, sim_fq_files[0] + '.gz'),
                context.write_intermediate_file(job, sim_fq_files[1] + '.gz')]
    else:
        sim_fq_file = os.path.join(work_dir, 'sim.fq.gz')
        cmd = [['samtools', 'fastq', os.path.basename(bam_file), '-N']]
        # we change /1 /2 --> _1 _2 to be compatible with rest of mapeval
        cmd.append(['sed', '-e', 's/\/1/_1/g', '-e', 's/\/2/_2/g'])
        cmd.append(['gzip'])
        with open(sim_fq_file, 'w') as sim_file:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = sim_file)
        return [context.write_intermediate_file(job, sim_fq_file)]
    
def run_gam_to_fastq(job, context, gam_file_id, paired_mode,
                     add_paired_suffix=False, out_name = 'sim', out_store = False, ):
    """
    convert a gam to fastq (or pair of fastqs)
    """
    work_dir = job.fileStore.getLocalTempDir()

    # read the gam file
    gam_file = os.path.join(work_dir, 'input.gam')
    job.fileStore.readGlobalFile(gam_file_id, gam_file)

    write_fn = context.write_output_file if out_store else context.write_intermediate_file
    
    # if we're paired, must make some split files
    if paired_mode:
        # convert to json (todo: have docker image that can do vg and jq)
        json_file = gam_file + '.json'
        cmd = ['vg', 'view', '-a', os.path.basename(gam_file)]
        with open(json_file, 'w') as out_json:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_json)

        sim_fq_files = [None, os.path.join(work_dir, '{}_1{}.fq.gz'.format(out_name, 's' if add_paired_suffix else '')),
                        os.path.join(work_dir, '{}_2{}.fq.gz'.format(out_name, 's' if add_paired_suffix else ''))]

        # make a fastq for each end of pair
        for i in [1, 2]:
            # extract paired end with jq
            cmd = ['jq', '-cr', 'select(.name | test("_{}$"))'.format(i),
                   os.path.basename(json_file)]
            end_file = json_file + '.{}'.format(i)
            with open(end_file, 'w') as end_out:
                context.runner.call(job, cmd, work_dir = work_dir, outfile = end_out)

            cmd = [['vg', 'view', '-JaG', os.path.basename(end_file)]]
            cmd.append(['vg', 'view', '-X', '-'])
            if not add_paired_suffix:
                cmd.append(['sed', 's/_{}$//'.format(i)])
            cmd.append(['gzip'])

            with open(sim_fq_files[i], 'w') as sim_out:
                context.runner.call(job, cmd, work_dir = work_dir, outfile = sim_out)

            os.remove(end_file)

        return [write_fn(job, sim_fq_files[1]), write_fn(job, sim_fq_files[2])]
            
    else:
        # extract reads from gam.  as above, need to have single docker container (which shouldn't be
        # a big deal) to run all these chained command and avoid huge files on disk
        extracted_reads_file = os.path.join(work_dir, '{}.fq.gz'.format(out_name))
        cmd = [['vg', 'view', '-X', os.path.basename(gam_file)]]
        cmd.append(['gzip'])
        with open(extracted_reads_file, 'w') as out_ext:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_ext)

        return [write_fn(job, extracted_reads_file)]

def run_concat_fastqs(job, context, fq_reads_ids):
    """ concatenate some fastq files
    """
    work_dir = job.fileStore.getLocalTempDir()

    assert len(fq_reads_ids) == 2
    # read the reads
    fq_file_names = [os.path.join(work_dir, 'reads-{}.fq.gz'.format(i)) \
                     for i in range(len(fq_reads_ids))]
    for fq_id, fq_name in zip(fq_reads_ids, fq_file_names):
        job.fileStore.readGlobalFile(fq_id, fq_name, mutable=fq_name==fq_file_names[0])

    # concat the reads (should work fine for gzipped or not)
    with open(fq_file_names[0], 'a') as out_file:
        for fq_name in fq_file_names[1:]:
            with open(fq_name) as fq_file:
                shutil.copyfileobj(fq_file, out_file)

    return context.write_intermediate_file(job, fq_file_names[0])

def run_strip_fq_ext(job, context, fq_reads_ids):
    """ bwa can't read reads with _1 _2 extensions for paired end alignment.  strip here
    """
    
    work_dir = job.fileStore.getLocalTempDir()

    # read the reads
    fq_file_names = [os.path.join(work_dir, 'reads-{}.fq.gz'.format(i)) \
                     for i in range(len(fq_reads_ids))]
    out_file_names = [os.path.join(work_dir, 'reads-strip-{}.fq.gz'.format(i)) \
                      for i in range(len(fq_reads_ids))]
    out_ids = []
    
    for fq_id, fq_name,  out_name in zip(fq_reads_ids, fq_file_names, out_file_names):
        job.fileStore.readGlobalFile(fq_id, fq_name, mutable=fq_name==fq_file_names[0])
        cmd = [['bgzip', '-dc', os.path.basename(fq_name)]]
        cmd.append(['sed', '-e', 's/_1$\|_2$//g'])
        cmd.append(['bgzip', '-c'])
        with open(out_name, 'w') as out_file:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_file)
        out_ids.append(context.write_intermediate_file(job, out_name))

    return out_ids
    
def run_bwa_mem(job, context, fq_reads_ids, bwa_index_ids, paired_mode):
    """ run bwa-mem on reads in a gam.  optionally run in paired mode
    return id of bam file
    """

    work_dir = job.fileStore.getLocalTempDir()

    # read the reads
    fq_file_names = []
    for i, fq_reads_id in enumerate(fq_reads_ids):
        fq_file_names.append(os.path.join(work_dir, 'reads{}.fq.gz'.format(i)))
        job.fileStore.readGlobalFile(fq_reads_id, fq_file_names[-1])

    # and the index files
    fasta_file = os.path.join(work_dir, 'reference.fa')
    for suf, idx_id in bwa_index_ids.items():
        job.fileStore.readGlobalFile(idx_id, '{}{}'.format(fasta_file, suf))

    # output positions file
    bam_file = os.path.join(work_dir, 'bwa-mem')
    if paired_mode:
        bam_file += '-pe'
    bam_file += '.bam'
    
    # if we're paired, must make some split files
    if paired_mode:

        # run bwa-mem on the paired end input
        start_time = timeit.default_timer()
        cmd = ['bwa', 'mem', '-t', str(context.config.alignment_cores), os.path.basename(fasta_file),
                os.path.basename(fq_file_names[0])]
        if len(fq_file_names) == 2:
            cmd += [os.path.basename(fq_file_names[1])]
        # if one file comes in, it had better be interleaved            
        else:
            cmd += ['-p']
        cmd += context.config.bwa_opts
        
        with open(bam_file + '.sam', 'w') as out_sam:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_sam)

        end_time = timeit.default_timer()
        run_time = end_time - start_time

        # we take care to mimic output message from vg_map.py, so we can mine them both for the jenkins
        # report        
        RealtimeLogger.info("Aligned aligned-linear_0.gam. Process took {} seconds with paired-end bwa-mem".format(
            run_time))
            
        # separate samtools for docker (todo find image with both)
        # 2304 = get rid of 256 (secondary) + 2048 (supplementary)        
        cmd = ['samtools', 'view', '-1', '-F', '2304', os.path.basename(bam_file + '.sam')]
        with open(bam_file, 'w') as out_bam:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_bam)

    # single end
    else:
        assert len(fq_file_names) == 1

        # run bwa-mem on single end input
        start_time = timeit.default_timer()
        cmd = ['bwa', 'mem', '-t', str(context.config.alignment_cores), os.path.basename(fasta_file),
                os.path.basename(fq_file_names[0])] + context.config.bwa_opts

        with open(bam_file + '.sam', 'w') as out_sam:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_sam)

        end_time = timeit.default_timer()
        run_time = end_time - start_time

        # we take care to mimic output message from vg_map.py, so we can mine them both for the jenkins
        # report
        RealtimeLogger.info("Aligned aligned-linear_0.gam. Process took {} seconds with single-end bwa-mem".format(
            run_time))            

        # separate samtools for docker (todo find image with both)
        # 2304 = get rid of 256 (secondary) + 2048 (supplementary)
        cmd = ['samtools', 'view', '-1', '-F', '2304', os.path.basename(bam_file + '.sam')]
        with open(bam_file, 'w') as out_bam:
            context.runner.call(job, cmd, work_dir = work_dir, outfile = out_bam)


    # return our id for the output bam file
    bam_file_id = context.write_output_file(job, bam_file)

    return bam_file_id, run_time

def extract_bam_read_stats(job, context, name, bam_file_id, paired, sep='_'):
    """
    extract positions, scores, and MAPQs from bam, return id of read stats file
    (lots of duplicated code with vg_sim, should merge?)
    
    Produces a read stats TSV of the format:
    read name, contig aligned to, alignment position, score, MAPQ
    
    TODO: Currently scores are not extracted and a score of 0 is always
    returned.

    """

    work_dir = job.fileStore.getLocalTempDir()

    # download input
    bam_file = os.path.join(work_dir, name)
    job.fileStore.readGlobalFile(bam_file_id, bam_file)

    out_pos_file = bam_file + '.tsv'

    # 2304 = get rid of 256 (secondary) + 2048 (supplementary)        
    cmd = [['samtools', 'view', os.path.basename(bam_file), '-F', '2304']]
    cmd.append(['grep', '-v', '^@'])
    if paired:
        # Now we use inline perl to parse the SAM flags and synthesize TSV
        # TODO: will need to switch to something more powerful to parse the score out of the AS tag. For now score everything as 0.
        # TODO: why _ and not / as the read name vs end number delimiter?
        # Note: we are now adding length/2 to the positions to be more consistent with vg annotate
        cmd.append(['perl', '-ne', '@val = split("\t", $_); print @val[0] . "{}" . (@val[1] & 64 ? "1" : @val[1] & 128 ? "2" : "?"), "\t" . @val[2] . "\t" . (@val[3] +  int(length(@val[9]) / 2)) . "\t0\t" . @val[4] . "\n";'.format(sep)])
    else:
        # No flags to parse since there's no end pairing and read names are correct.
        # Use inline perl again and insert a fake 0 score column
        # Note: we are now adding length/2 to the positions to be more consistent with vg annotate        
        cmd.append(['perl', '-ne', '@val = split("\t", $_); print @val[0] . "\t" . @val[2] . "\t" . (@val[3] +  int(length(@val[9]) / 2)) . "\t0\t" . @val[4] . "\n";'])
    cmd.append(['sort'])
    
    with open(out_pos_file, 'w') as out_pos:
        context.runner.call(job, cmd, work_dir = work_dir, outfile = out_pos)

    stats_file_id = context.write_intermediate_file(job, out_pos_file)
    return stats_file_id

    
def annotate_gam(job, context, xg_file_id, gam_file_id):
    """
    Annotate the given GAM file with positions from the given XG file.
    """
    
    work_dir = job.fileStore.getLocalTempDir()

    # download input
    RealtimeLogger.info('Download XG from file {}'.format(xg_file_id))
    xg_file = os.path.join(work_dir, 'index.xg')
    job.fileStore.readGlobalFile(xg_file_id, xg_file)
    gam_file = os.path.join(work_dir, 'reads.gam')
    job.fileStore.readGlobalFile(gam_file_id, gam_file)
    
    annotated_gam_file = os.path.join(work_dir, 'annotated.gam')
    
    cmd = [['vg', 'annotate', '-p', '-a', os.path.basename(gam_file), '-x', os.path.basename(xg_file)]]
    with open(annotated_gam_file, 'w') as out_file:
        context.runner.call(job, cmd, work_dir=work_dir, outfile=out_file)
    
    return context.write_intermediate_file(job, annotated_gam_file)
    
    
def extract_gam_read_stats(job, context, name, gam_file_id):
    """
    extract positions, scores, and MAPQs for reads from a gam, and return the id
    of the resulting read stats file
    
    Produces a read stats TSV of the format:
    read name, contig aligned to, alignment position, score, MAPQ
    
    If the GAM is not annotated with alignment positions, contig and position
    will both contain only "0" values.

    """

    work_dir = job.fileStore.getLocalTempDir()

    gam_file = os.path.join(work_dir, name)
    job.fileStore.readGlobalFile(gam_file_id, gam_file)

    out_pos_file = gam_file + '.tsv'
                           
    # go through intermediate json file until docker worked out
    gam_annot_json = gam_file + '.json'
    cmd = [['vg', 'view', '-aj', os.path.basename(gam_file)]]
    with open(gam_annot_json, 'w') as output_annot_json:
        context.runner.call(job, cmd, work_dir = work_dir, outfile=output_annot_json)

    # turn the annotated gam json into truth positions, as separate command since
    # we're going to use a different docker container.  (Note, would be nice to
    # avoid writing the json to disk)        
    jq_cmd = [['jq', '-c', '-r', '[.name, '
               'if .refpos != null then (.refpos[] | .name, .offset) else (null, null) end, '
               '.score, '
               'if .mapping_quality == null then 0 else .mapping_quality end ] | @tsv',
               os.path.basename(gam_annot_json)]]
    # convert back to _1 format (only relevant if running on bam input reads where / added automatically)
    jq_cmd.append(['sed', '-e', 's/null/0/g',  '-e', 's/\/1/_1/g', '-e', 's/\/2/_2/g'])

    with open(out_pos_file + '.unsorted', 'w') as out_pos:
        context.runner.call(job, jq_cmd, work_dir = work_dir, outfile=out_pos)

    # get rid of that big json asap
    os.remove(gam_annot_json)

    # sort the read stats file (not piping due to memory fears)
    sort_cmd = ['sort', os.path.basename(out_pos_file) + '.unsorted']
    with open(out_pos_file, 'w') as out_pos:
        context.runner.call(job, sort_cmd, work_dir = work_dir, outfile = out_pos)

    # Make sure each line has all columns
    RealtimeLogger.info("Make sure all lines are full length")
    context.runner.call(job, ['awk', '!length($5)',  os.path.basename(out_pos_file)], work_dir = work_dir)

    out_stats_file_id = context.write_intermediate_file(job, out_pos_file)
    return out_stats_file_id
    
def compare_positions(job, context, truth_file_id, name, stats_file_id, mapeval_threshold):
    """
    this is essentially pos_compare.py from vg/scripts
    return output file id.
    
    Compares positions from two TSV files. The truth has the format:
    read name, contig simulated from, true position
    
    And the file under test is a read stats TSV with the format:
    read name, contig aligned to, alignment position, alignment score, MAPQ
    
    Produces a CSV (NOT TSV) of the form:
    read name, correctness flag (0/1), MAPQ
    
    mapeval_threshold is the distance within which a mapping is held to have hit
    the correct position.
    
    """
    work_dir = job.fileStore.getLocalTempDir()

    true_read_stats_file = os.path.join(work_dir, 'true.tsv')
    job.fileStore.readGlobalFile(truth_file_id, true_read_stats_file)
    test_read_stats_file = os.path.join(work_dir, name + '.tsv')
    job.fileStore.readGlobalFile(stats_file_id, test_read_stats_file)

    out_file = os.path.join(work_dir, name + '.compare.positions')

    with open(true_read_stats_file) as truth, open(test_read_stats_file) as test, \
         open(out_file, 'w') as out:
        line_no = 0
        for true_fields, test_fields in itertools.izip(tsv.TsvReader(truth), tsv.TsvReader(test)):
            # Zip everything up and assume that the reads correspond
            line_no += 1
            # every input has a true position
            true_read_name = true_fields[0]
            if len(true_fields) < 3:
                # This seems to come up about one-in-a-million times from vg annotate as called
                # by toil-vg sim.  Once it is fixed, we can turn this back into an error
                logger.warning('Incorrect (< 3) true field counts on line {} for {}: {} and {}'.format(
                    line_no, name, true_fields, test_fields))
                true_fields = [true_read_name, '0', '0']
            if len(test_fields) < 5:
                # With the new TSV reader, the files should always have the
                # correct field counts. Some fields just might be empty.
                raise RuntimeError('Incorrect (<5) test field counts on line {} for {}: {} and {}'.format(
                    line_no, name, true_fields, test_fields))

            # map seq name->position
            true_pos_dict = dict(zip(true_fields[1::2], map(parse_int, true_fields[2::2])))
            aln_read_name = test_fields[0]
            if aln_read_name !=true_read_name:
                raise RuntimeError('Mismatch on line {} of {} and {}.  Read names differ: {} != {}'.format(
                    line_no, true_read_stats_file, test_read_stats_file, true_read_name, aln_read_name))
            aln_pos_dict = dict(zip(test_fields[1:-2:2], map(parse_int, test_fields[2:-2:2])))
            # Skip over score field
            aln_mapq = parse_int(test_fields[-1])
            aln_correct = 0
            for aln_chr, aln_pos in aln_pos_dict.items():
                if aln_chr in true_pos_dict and abs(true_pos_dict[aln_chr] - aln_pos) < mapeval_threshold:
                    aln_correct = 1
                    break

            out.write('{}, {}, {}\n'.format(aln_read_name, aln_correct, aln_mapq))
        
        # make sure same length
        has_next = False
        try:
            iter(truth).next()
            has_next = True
        except:
            pass
        try:
            iter(test).next()
            has_next = True
        except:
            pass
        if has_next:
            raise RuntimeError('read stats files have different lengths')
        
    out_file_id = context.write_intermediate_file(job, out_file)
    return out_file_id
    
def compare_scores(job, context, baseline_name, baseline_file_id, name, score_file_id):
    """
    Compares scores from TSV files. The baseline and file under test both have
    the format:
    read name, contig aligned to, alignment position, alignment score, MAPQ
    
    Produces a CSV (NOT TSV) of the form:
    read name, score difference, aligned score, baseline score
    
    If saved to the out store it will be:
    <condition name>.compare.<baseline name>.scores
    
    Uses the given (condition) name as a file base name for the file under test.
    
    """
    work_dir = job.fileStore.getLocalTempDir()

    baseline_read_stats_file = os.path.join(work_dir, 'baseline.tsv')
    job.fileStore.readGlobalFile(baseline_file_id, baseline_read_stats_file)
    test_read_stats_file = os.path.join(work_dir, name + '.tsv')
    job.fileStore.readGlobalFile(score_file_id, test_read_stats_file)

    out_file = os.path.join(work_dir, '{}.compare.{}.scores'.format(name, baseline_name))

    with open(baseline_read_stats_file) as baseline, open(test_read_stats_file) as test:
        with open(out_file, 'w') as out:
            line_no = 0
            for baseline_fields, test_fields in itertools.izip(tsv.TsvReader(baseline), tsv.TsvReader(test)):
                # Zip everything up and assume that the reads correspond
                line_no += 1
                
                if len(baseline_fields) < 5 or len(test_fields) < 5:
                    raise RuntimeError('Incorrect field counts on line {} for {}: {} and {}'.format(
                        line_no, name, baseline_fields, test_fields))
                
                if baseline_fields[0] != test_fields[0]:
                    # Read names must correspond or something has gone wrong
                    raise RuntimeError('Mismatch on line {} of {} and {}.  Read names differ: {} != {}'.format(
                        line_no, baseline_read_stats_file, test_read_stats_file, baseline_fields[0], test_fields[0]))
                
                # Order is: name, conting, pos, score, mapq
                aligned_score = test_fields[-2]
                baseline_score = baseline_fields[-2]
                # Compute the score difference. Scores are integers.
                score_diff = parse_int(aligned_score) - parse_int(baseline_score)
                
                # Report the score difference            
                out.write('{}, {}, {}, {}\n'.format(baseline_fields[0], score_diff, aligned_score, baseline_score))
        
        # Save stats file for inspection
        out_file_id = context.write_intermediate_file(job, out_file)
        
        # make sure same length
        has_next = False
        found_line1 = None
        found_line2 = None
        try:
            found_line1 = iter(baseline).next().rstrip()
            has_next = True
        except:
            pass
        try:
            found_line2 = iter(test).next().rstrip()
            has_next = True
        except:
            pass
        if has_next:
            raise RuntimeError('read stats files have different lengths ({}, {})'.format(found_line1, found_line2))
        
    return out_file_id

def run_map_eval_index(job, context, xg_file_ids, gcsa_file_ids, gbwt_file_ids, id_range_file_ids, snarl_file_ids, vg_file_ids):
    """ 
    Index the given vg files.
    
    If no vg files are provided, pass through the given indexes. Indexes are
    lists of index IDs, one per graph, except gcsa_file_ids, which is tuples of
    GCSA and LCP file IDs, one tuple per graph. Index types which are not used
    should have falsey values instead of lists.
    
    Returns a list of dicts from index type name to index file ID, as used by
    run_indexing in vg_index.py, holding file IDs for different index
    components.
    
    """

    # index_ids are dicts from index type to file ID as returned by run_indexing
    index_ids = []
    if vg_file_ids:
        for vg_file_id in vg_file_ids:
            index_job = job.addChildJobFn(run_indexing, context, [vg_file_id], ['default.vg'],
                                          'index', ['default'], 
                                          cores=context.config.misc_cores, memory=context.config.misc_mem,
                                          disk=context.config.misc_disk)
            index_ids.append(index_job.rv())
    else:
        for i, xg_id in enumerate(xg_file_ids):
            # For each graph, gather and tag its indexes
            indexes = {}
            indexes['xg'] = xg_id
            if gcsa_file_ids:
                indexes['gcsa'], indexes['lcp'] = gcsa_file_ids[i]
            if gbwt_file_ids:
                indexes['gbwt'] = gbwt_file_ids[i]
            if id_range_file_ids:
                indexes['id_ranges'] = id_range_file_ids[i]
            if snarl_file_ids:
                indexes['snarls'] = snarl_file_ids[i]
                
            # Put the indexes in the list of index dicts for each graph
            index_ids.append(indexes)

    
    return index_ids

def run_map_eval_align(job, context, index_ids, gam_names, gam_file_ids,
                       reads_fastq_single_ids, reads_fastq_paired_ids, reads_fastq_paired_for_vg_ids,
                       fasta_file_id, bwa_index_ids, do_bwa, do_single,
                       do_paired, singlepath, multipath, ignore_quals):
    """
    Run alignments, if alignment files have not already been provided.
    
    Returns a list of graph/gam names, a list of associated gam file IDs, a
    list of associated xg index IDs, and a list of BAM file IDs (or Nones) for
    realigned read BAMs, and a list of running times (for map commands not including toil-vg overhead)
    
    We need to modify the name and index lists because we synthesize paired-end
    versions of existing entries.
    
    """

    # scrape out the xg ids, don't need others any more after this step
    xg_ids = [index_id['xg'] for index_id in index_ids]

    # the ids and names we pass forward
    out_xg_ids = xg_ids if gam_file_ids else []
    out_gam_names = gam_names if gam_file_ids else []

    # the map times
    map_times = [None] * len(gam_file_ids) if gam_file_ids else []

    # we use this hack to run multiple batches of mpmap opts
    mpmap_opts_list = [context.config.mpmap_opts]
    if context.config.more_mpmap_opts:
        mpmap_opts_list += context.config.more_mpmap_opts

    if ignore_quals:
        # Make sure we don't use quality adjusted alignment since simulation doesn't make qualities
        for mpmap_opts in mpmap_opts_list:
            if '-A' not in mpmap_opts and '--no-qual-adjust' not in mpmap_opts:
                mpmap_opts.append('-A')
        context.config.map_opts = [o for o in context.config.map_opts if o not in ['-A', '--qual-adjust']]

    assert reads_fastq_single_ids or reads_fastq_paired_ids
    def fq_names(fq_reads_ids):
        return ['input{}.fq.gz'.format(i) for i in range(len(fq_reads_ids))]

    do_vg_mapping = not gam_file_ids and (singlepath or multipath)
    if do_vg_mapping and singlepath and do_single:
        gam_file_ids = []
        # run vg map if requested
        for i, indexes in enumerate(index_ids):
            map_job = job.addChildJobFn(run_mapping, context, fq_names(reads_fastq_single_ids),
                                        None, None, 'aligned-{}'.format(gam_names[i]),
                                        False, False, indexes,
                                        reads_fastq_single_ids,
                                        cores=context.config.misc_cores,
                                        memory=context.config.misc_mem, disk=context.config.misc_disk)
            gam_file_ids.append(map_job.rv(0))
            map_times.append(map_job.rv(1))

        # make sure associated lists are extended
        out_xg_ids += xg_ids
        out_gam_names += gam_names

    # Do the single-ended multipath mapping
    if do_vg_mapping and multipath and do_single:
        for opt_num, mpmap_opts in enumerate(mpmap_opts_list):
            mp_context = copy.deepcopy(context)
            mp_context.config.mpmap_opts = mpmap_opts
            for i, indexes in enumerate(index_ids):
                map_job = job.addChildJobFn(run_mapping, mp_context, fq_names(reads_fastq_single_ids),
                                            None, None, 'aligned-{}-mp'.format(gam_names[i]),
                                            False, True, indexes,
                                            reads_fastq_single_ids,
                                            cores=mp_context.config.misc_cores,
                                            memory=mp_context.config.misc_mem, disk=mp_context.config.misc_disk)
                gam_file_ids.append(map_job.rv(0))
                map_times.append(map_job.rv(1))

            # make sure associated lists are extended to fit new paired end mappings
            out_xg_ids += xg_ids
            out_gam_names += [n + '-mp{}'.format(opt_num if opt_num > 0 else '') for n in gam_names]

    if do_vg_mapping and do_paired and singlepath:
        # run paired end version of all vg inputs if --pe-gams specified
        for i, indexes in enumerate(index_ids):
            interleaved = len(reads_fastq_paired_for_vg_ids) == 1
            map_job = job.addChildJobFn(run_mapping, context, fq_names(reads_fastq_paired_for_vg_ids),
                                        None, None, 'aligned-{}-pe'.format(gam_names[i]),
                                        interleaved, False, indexes,
                                        reads_fastq_paired_for_vg_ids,
                                        cores=context.config.misc_cores,
                                        memory=context.config.misc_mem, disk=context.config.misc_disk)
            gam_file_ids.append(map_job.rv(0))
            map_times.append(map_job.rv(1))            
            
        # make sure associated lists are extended to fit new paired end mappings
        out_xg_ids += xg_ids
        out_gam_names += [n + '-pe' for n in gam_names]

    # Do the paired-ended multipath mapping
    if do_vg_mapping and do_paired and multipath:
        interleaved = len(reads_fastq_paired_for_vg_ids) == 1
        for opt_num, mpmap_opts in enumerate(mpmap_opts_list):
            mp_context = copy.deepcopy(context)
            mp_context.config.mpmap_opts = mpmap_opts
            for i, indexes in enumerate(index_ids):
                map_job = job.addChildJobFn(run_mapping, mp_context, fq_names(reads_fastq_paired_for_vg_ids),
                                            None, None, 'aligned-{}-mp-pe'.format(gam_names[i]),
                                            interleaved, True, indexes,
                                            reads_fastq_paired_for_vg_ids,
                                            cores=mp_context.config.misc_cores,
                                            memory=mp_context.config.misc_mem, disk=mp_context.config.misc_disk)
                gam_file_ids.append(map_job.rv(0))
                map_times.append(map_job.rv(1))            

            # make sure associated lists are extended to fit new paired end mappings
            out_xg_ids += xg_ids
            out_gam_names += [n + '-mp{}-pe'.format(opt_num if opt_num > 0 else '') for n in gam_names]
    
    # run bwa if requested
    bwa_bam_file_ids, bwa_mem_times = [None, None], [None, None]
    if do_bwa:
        bwa_start_job = Job()
        job.addChild(bwa_start_job)
        bwa_index_job = bwa_start_job.addChildJobFn(run_bwa_index, context,
                                                    fasta_file_id, bwa_index_ids,
                                                    cores=context.config.alignment_cores, memory=context.config.alignment_mem,
                                                    disk=context.config.alignment_disk)
        bwa_index_ids = bwa_index_job.rv()
                
        if do_single:
            bwa_mem_job = bwa_start_job.addFollowOnJobFn(run_bwa_mem, context, reads_fastq_single_ids, bwa_index_ids, False,
                                                         cores=context.config.alignment_cores, memory=context.config.alignment_mem,
                                                         disk=context.config.alignment_disk)
            bwa_bam_file_ids[0] = bwa_mem_job.rv(0)
            bwa_mem_times[0] = bwa_mem_job.rv(1)
        if do_paired:
            bwa_mem_job = bwa_start_job.addFollowOnJobFn(run_bwa_mem, context, reads_fastq_paired_ids, bwa_index_ids, True,
                                                         cores=context.config.alignment_cores, memory=context.config.alignment_mem,
                                                         disk=context.config.alignment_disk)
            bwa_bam_file_ids[1] = bwa_mem_job.rv(0)
            bwa_mem_times[1] = bwa_mem_job.rv(1)

    return out_gam_names, gam_file_ids, out_xg_ids, map_times, bwa_bam_file_ids, bwa_mem_times
    
def run_map_eval_comparison(job, context, xg_file_ids, gam_names, gam_file_ids,
                            bam_names, bam_file_ids, pe_bam_names, pe_bam_file_ids,
                            bwa_bam_file_ids, true_read_stats_file_id, mapeval_threshold,
                            score_baseline_name=None, original_read_gam=None):
    """
    run the mapping comparison.  Dump some tables into the outstore.
    
    Returns a pair of the position comparison results and the score comparison results.
    
    The score comparison results are a dict from baseline name to comparison
    against that baseline. Each comparison's data is a tuple of a list of
    individual per-graph comparison file IDs and an overall stats file for that
    comparison.
    
    If score_baseline_name is specified, all GAMs have their scores compared
    against the scores from the GAM with that name as a baseline.
    
    If original_read_gam is specified, all GAMs have their scores compared
    against that GAM's scores as a baseline.
    
    Each result set is itself a pair, consisting of a list of per-graph file IDs, and an overall statistics file ID.
    
    """
    
    # munge out the returned pair from run_bwa_index()
    if bwa_bam_file_ids[0] is not None:
        bam_file_ids.append(bwa_bam_file_ids[0])
        bam_names.append('bwa-mem')
    if bwa_bam_file_ids[1] is not None:
        pe_bam_file_ids.append(bwa_bam_file_ids[1])
        pe_bam_names.append('bwa-mem-pe')

    # We need to keep the BAM stats jobs around to wait on them
    bam_stats_jobs = []

    # get the bwa read alignment statistics, one id for each bam_name
    bam_stats_file_ids = []
    for bam_i, bam_id in enumerate(bam_file_ids):
        name = '{}-{}.bam'.format(bam_names[bam_i], bam_i)
        bam_stats_jobs.append(job.addChildJobFn(extract_bam_read_stats, context, name, bam_id, False,
                                                cores=context.config.misc_cores, memory=context.config.misc_mem,
                                                disk=context.config.misc_disk))
        bam_stats_file_ids.append(bam_stats_jobs[-1].rv())
    # separate flow for paired end bams because different logic used
    pe_bam_stats_file_ids = []
    for bam_i, bam_id in enumerate(pe_bam_file_ids):
        name = '{}-{}.bam'.format(pe_bam_names[bam_i], bam_i)
        bam_stats_jobs.append(job.addChildJobFn(extract_bam_read_stats, context, name, bam_id, True,
                                                cores=context.config.misc_cores, memory=context.config.misc_mem,
                                                disk=context.config.misc_disk))
        pe_bam_stats_file_ids.append(bam_stats_jobs[-1].rv())

    # get the gam read alignment statistics, one for each gam_name (todo: run vg map like we do bwa?)
    gam_stats_file_ids = []
    # We also need to keep the jobs around so we can enforce that things run after them.
    gam_stats_jobs = []
    for gam_i, gam_id in enumerate(gam_file_ids):
        name = '{}-{}.gam'.format(gam_names[gam_i], gam_i)
        # run_mapping will return a list of gam_ids.  since we don't
        # specify id ranges, this will always have one entry
        gam = gam_id
        if type(gam_id) is list:
            assert len(gam_id) == 1
            gam = gam_id[0]
        RealtimeLogger.info('Work on GAM {} = {} named {} with XG id {}'.format(gam_i, gam, gam_names[gam_i], xg_file_ids[gam_i]))
        
        # Make a job to annotate the GAM
        annotate_job = job.addChildJobFn(annotate_gam, context, xg_file_ids[gam_i], gam,
                                         cores=context.config.misc_cores, memory=context.config.misc_mem,
                                         disk=context.config.misc_disk)
        
        # Then compute stats on the annotated GAM
        gam_stats_jobs.append(annotate_job.addFollowOnJobFn(extract_gam_read_stats, context,
                                                            name, annotate_job.rv(),
                                                            cores=context.config.misc_cores, memory=context.config.misc_mem,
                                                            disk=context.config.misc_disk))
        
        gam_stats_file_ids.append(gam_stats_jobs[-1].rv())

    # compare all our positions, and dump results to the out store. Get a tuple
    # of individual comparison files and overall stats file.
    position_comparison_job = job.addChildJobFn(run_map_eval_compare_positions, context,
                                                true_read_stats_file_id, gam_names, gam_stats_file_ids,
                                                bam_names, bam_stats_file_ids, pe_bam_names, pe_bam_stats_file_ids,
                                                mapeval_threshold,
                                                cores=context.config.misc_cores, memory=context.config.misc_mem,
                                                disk=context.config.misc_disk)
    for dependency in itertools.chain(gam_stats_jobs, bam_stats_jobs):
        dependency.addFollowOn(position_comparison_job)
    position_comparison_results = position_comparison_job.rv()
    
    # This will map from baseline name to score comparison data against that
    # baseline
    score_comparisons = {}
    
    if score_baseline_name is not None:
        # We want to compare the scores from all the GAMs against a baseline
        
        # Make a dict mapping from assigned GAM name in gam_names to the stats file for that GAM's alignment
        name_to_stats_id = dict(itertools.izip(gam_names, gam_stats_file_ids))
        
        # Find the baseline scores
        baseline_stats_id = name_to_stats_id[score_baseline_name]
        
        # compare all our scores against the baseline, and dump results to the
        # out store. 
        score_comp_job = job.addChildJobFn(run_map_eval_compare_scores, context, score_baseline_name, baseline_stats_id,
                                           gam_names, gam_stats_file_ids, bam_names, bam_stats_file_ids,
                                           pe_bam_names, pe_bam_stats_file_ids, cores=context.config.misc_cores,
                                           memory=context.config.misc_mem, disk=context.config.misc_disk)
        for dependency in itertools.chain(gam_stats_jobs, bam_stats_jobs):
            dependency.addFollowOn(score_comp_job)
                             
        # Get a tuple of individual comparison files and overall stats file.
        score_comparisons[score_baseline_name] = score_comp_job.rv()
        
    if original_read_gam is not None:
        # Also compare against the original GAM's scores as a baseline
        
        # First compute its stats file
        stats_job = job.addChildJobFn(extract_gam_read_stats, context,
                                      name, original_read_gam,
                                      cores=context.config.misc_cores, memory=context.config.misc_mem,
                                      disk=context.config.misc_disk)
        
        # compare all our scores against this other baseline, and dump results
        # to the out store.
        score_comp_job = stats_job.addFollowOnJobFn(run_map_eval_compare_scores, context, 'input', stats_job.rv(),
                                                    gam_names, gam_stats_file_ids, bam_names, bam_stats_file_ids,
                                                    pe_bam_names, pe_bam_stats_file_ids, cores=context.config.misc_cores,
                                                    memory=context.config.misc_mem, disk=context.config.misc_disk)
                                                    
        for dependency in itertools.chain(gam_stats_jobs, bam_stats_jobs):
            dependency.addFollowOn(score_comp_job)
                                                    
        # Save the results
        score_comparisons['input'] = score_comp_job.rv()
        
        
        
    return position_comparison_results, score_comparisons

def run_map_eval_compare_positions(job, context, true_read_stats_file_id, gam_names, gam_stats_file_ids,
                         bam_names, bam_stats_file_ids, pe_bam_names, pe_bam_stats_file_ids, mapeval_threshold):
    """
    Compare the read positions for each read across the different aligmment
    methods.
    
    Produces a bunch of individual comparison files against the truth, a
    combined "positions.results.tsv" across all aligners, and a statistics file
    "stats.tsv" in the out_store.
    
    Returns the list of comparison files, and the stats file ID.
    """

    # merge up all the output data into one list
    names = gam_names + bam_names + pe_bam_names
    stats_file_ids = gam_stats_file_ids + bam_stats_file_ids + pe_bam_stats_file_ids

    compare_ids = []
    for name, stats_file_id in zip(names, stats_file_ids):
        compare_ids.append(job.addChildJobFn(compare_positions, context, true_read_stats_file_id, name,
                                             stats_file_id, mapeval_threshold,
                                             cores=context.config.misc_cores, memory=context.config.misc_mem,
                                             disk=context.config.misc_disk).rv())

    position_comp_file_id = job.addFollowOnJobFn(run_process_position_comparisons, context, names, compare_ids,
                                            cores=context.config.misc_cores, memory=context.config.misc_mem,
                                            disk=context.config.misc_disk).rv(1)
                                         
    return compare_ids, position_comp_file_id

def run_process_position_comparisons(job, context, names, compare_ids):
    """
    Write some raw tables of position comparisons to the output.  Compute some
    stats for each graph.
    
    Returns (the stats file's file ID, the position results file's ID)
    """

    work_dir = job.fileStore.getLocalTempDir()

    map_stats = []

    # make the position.results.tsv and position.stats.tsv
    results_file = os.path.join(work_dir, 'position.results.tsv')
    with open(results_file, 'w') as out_results:
        out_results.write('correct\tmq\taligner\tread\n')

        def write_tsv(comp_file, a):
            """
            Read the given comparison CSV for the given condition name, and dump
            it to the combined results file.
            """
            # tack on '-se' to single endings to make more consistent ordering in r plots
            method = a + '-se' if not (a.endswith('-pe') or a.endswith('-se')) else a
            
            with open(comp_file) as comp_in:
                for line in comp_in:
                    toks = line.rstrip().split(', ')
                    out_results.write('{}\t{}\t{}\t{}\n'.format(toks[1], toks[2], method, toks[0]))

        for name, compare_id in itertools.izip(names, compare_ids):
            compare_file = os.path.join(work_dir, '{}.compare.positions'.format(name))
            job.fileStore.readGlobalFile(compare_id, compare_file)
            context.write_output_file(job, compare_file)
            write_tsv(compare_file, name)

            map_stats.append([job.addChildJobFn(run_acc, context, name, compare_id, cores=context.config.misc_cores,
                                                memory=context.config.misc_mem, disk=context.config.misc_disk).rv(),
                              job.addChildJobFn(run_auc, context, name, compare_id, cores=context.config.misc_cores,
                                                memory=context.config.misc_mem, disk=context.config.misc_disk).rv(),
                              job.addChildJobFn(run_qq, context, name, compare_id, cores=context.config.misc_cores,
                                                memory=context.config.misc_mem, disk=context.config.misc_disk).rv(),
                              job.addChildJobFn(run_max_f1, context, name, compare_id, cores=context.config.misc_cores,
                                                memory=context.config.misc_mem, disk=context.config.misc_disk).rv()])
            
    position_results_id = context.write_output_file(job, results_file)

    return job.addFollowOnJobFn(run_write_position_stats, context, names, map_stats).rv(), position_results_id

def run_write_position_stats(job, context, names, map_stats):
    """
    write the position comparison statistics as tsv, both to the Toil fileStore
    and to the out_store as "stats.tsv".
    
    Returns the ID of the file written.
    
    This is different than the stats TSV format used internally, for read stats.
    """

    work_dir = job.fileStore.getLocalTempDir()
    stats_file = os.path.join(work_dir, 'stats.tsv')
    with open(stats_file, 'w') as stats_out:
        stats_out.write('aligner\tcount\tacc\tauc\tqq-r\tmax-f1\n')
        for name, stats in zip(names, map_stats):
            stats_out.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(name, stats[0][0], stats[0][1],
                                                          stats[1][0], stats[2], stats[3]))

    stats_file_id = context.write_output_file(job, stats_file)
    
    return stats_file_id
    
def run_acc(job, context, name, compare_id):
    """
    Percentage of correctly aligned reads (ignore quality)
    """
    
    work_dir = job.fileStore.getLocalTempDir()

    compare_file = os.path.join(work_dir, '{}.compare.positions'.format(name))
    job.fileStore.readGlobalFile(compare_id, compare_file)
    
    total = 0
    correct = 0
    with open(compare_file) as compare_f:
        for line in compare_f:
            toks = line.split(', ')
            total += 1
            if toks[1] == '1':
                correct += 1
                
    acc = float(correct) / float(total) if total > 0 else 0
    return total, acc
    
def run_auc(job, context, name, compare_id):
    """
    AUC of roc plot.
    
    ROC plot is defined with mismapped reads being negatives, correctly-mapped
    reads being positives, and AUC expressing how good of a classifier of
    correctly-mapped-ness the MAPQ score is. It says nothing about how well the
    reads are actually mapped.
    
    """
    if not have_sklearn:
        return ["sklearn_not_installed"] * 2 
    
    work_dir = job.fileStore.getLocalTempDir()

    compare_file = os.path.join(work_dir, '{}.compare.positions'.format(name))
    job.fileStore.readGlobalFile(compare_id, compare_file)

    try:
        data = np.loadtxt(compare_file, dtype=np.int, delimiter =', ', usecols=(1,2)).T
        auc = roc_auc_score(data[0], data[1])
        aupr = average_precision_score(data[0], data[1])
    except:
        # will happen if file is empty
        auc, aupr = 0, 0

    return auc, aupr
    
def run_max_f1(job, context, name, compare_id):
    """
    Compute and return maximum F1 score for correctly mapping reads, using MAPQ as a confidence.
    
    The problem is that read mapping is a multi-class classification problem with ~ 1 class per genome base, and F1 is a 2-class classification statistic. So we squish the concept a bit.
    
    For a given MAPQ threshold, we define a fake confusion matrix:
    
    TP = reads meeting threshold mapped correctly
    FP = reads meeting threshold mapped incorrectly
    TN = reads that didn't come from the target graph and didn't meet the threshold (i.e. 0 reads)
    FN = reads that weren't both correctly mapped and assigned a MAPQ meeting the threshold
    
    Then we calculate precision = TP / (TP + FP) and recall = TP / (TP + FN), and from those calculate an F1.
    Then we calculate the best F1 across all the MAPQ values.
    
    """
    if not have_sklearn:
        return ["sklearn_not_installed"] * 2 
    
    work_dir = job.fileStore.getLocalTempDir()

    compare_file = os.path.join(work_dir, '{}.compare.positions'.format(name))
    job.fileStore.readGlobalFile(compare_id, compare_file)

    # Load up the correct/incorrect flag (data[_, 0]) and the scores (data[_, 1])
    data = np.loadtxt(compare_file, dtype=np.int, delimiter =', ', usecols=(1,2))
    
    # Sort on score (see <https://stackoverflow.com/a/2828121/402891>) in
    # descending order. So reads we want to take first come first.
    data = data[data[:,1].argsort()[::-1]]
    
    # What's the last MAPQ we did?
    last_mapq = None
    
    # What so far is our confusion matrix?
    tp = 0
    fp = 0
    tn = 0
    fn = len(data)
    
    # What's our max f1?
    max_f1 = 0
    
    # And how do we calculate it?
    def emit_f1():
        if tp > 0 or (fp > 0 and fn > 0):
            # Safe to compute precision and recall, so do that.
            precision = float(tp) / (tp + fp)
            recall = float(tp) / (tp + fn)
            if precision > 0 or recall > 0:
                # Safe to compute an F1, so do that
                f1 = 2 * (precision * recall) / (precision + recall)
                return max(max_f1, f1)
        return max_f1
    
    for correct, score in data:
        # For each read in descending MAPQ order
        
        if score != last_mapq:
            # We're moving on to a new tranche
            max_f1 = emit_f1()
            
        # This read is now a positive. It may be true or false.
        if correct:
            tp += 1
        else:
            fp += 1
        # It is no longer a false negative
        fn -= 1
        
    # At the end of all the reads, that's another tranche done
    max_f1 = emit_f1()
        
    return max_f1

def run_qq(job, context, name, compare_id):
    """
    some measure of qq consistency 
    """
    if not have_sklearn:
        return "sklearn_not_installed"

    work_dir = job.fileStore.getLocalTempDir()

    compare_file = os.path.join(work_dir, '{}.compare.positions'.format(name))
    job.fileStore.readGlobalFile(compare_id, compare_file)

    try:
        data = np.loadtxt(compare_file, dtype=np.int, delimiter =', ', usecols=(1,2))

        # this can surley be sped up if necessary
        correct = Counter()
        total = Counter()
        for row in data:
            correct[row[1]] += row[0]
            total[row[1]] += 1

        qual_scores = []
        qual_observed = []            
        for qual, cor in correct.items():
            qual_scores.append(qual)
            p_err = max(1. - float(cor) / float(total[qual]), sys.float_info.epsilon)
            observed_score =-10. * math.log10(p_err)
            qual_observed.append(observed_score)

        # should do non-linear regression as well? 
        r2 = r2_score(qual_observed, qual_scores)
    except:
        # will happen if file is empty
        r2 = 'fail'

    return r2
    
def run_map_eval_compare_scores(job, context, baseline_name, baseline_stats_file_id, gam_names, gam_stats_file_ids,
                                bam_names, bam_stats_file_ids, pe_bam_names, pe_bam_stats_file_ids):
    """
    Compare scores in the given stats files in the lists to those in the given
    baseline stats file.
    
    Stats file format is a TSV of:
    read name, contig name, contig offset, score, mapping quality
    
    Will save the output to the outstore, as a CSV of read name and score
    difference, named <GAM/BAM name>.compare.<baseline name>.scores.
    
    Will also save a concatenated TSV file, with score difference and quoted
    aligner/condition name, as score.results.<baseline_name>.tsv
    
    Returns the list of comparison file IDs and the score results file ID.
    
    For now, just ignores BAMs because we don't pull in pysam to parse out their
    scores.
    
    """
    
    # merge up all the condition names and file IDs into synchronized lists
    # TODO: until we can extract the scores from BAMs, just process the GAMs
    names = gam_names
    stats_file_ids = gam_stats_file_ids
    
    RealtimeLogger.info(names)
    RealtimeLogger.info(stats_file_ids)

    compare_ids = []
    for name, stats_file_id in zip(names, stats_file_ids):
        compare_ids.append(job.addChildJobFn(compare_scores, context, baseline_name, baseline_stats_file_id, name, stats_file_id,
                                             cores=context.config.misc_cores, memory=context.config.misc_mem,
                                             disk=context.config.misc_disk).rv())

    stats_job = job.addFollowOnJobFn(run_process_score_comparisons, context, baseline_name, names, compare_ids,
                                     cores=context.config.misc_cores, memory=context.config.misc_mem,
                                     disk=context.config.misc_disk)
                                     
    return compare_ids, stats_job.rv()

def run_process_score_comparisons(job, context, baseline_name, names, compare_ids):
    """
    Write some raw tables of score comparisons against the given baseline to the
    output.  Compute some stats for each graph.
    
    Returns the file ID of the overall stats file "score.stats.<baseline name>.tsv".
    """

    work_dir = job.fileStore.getLocalTempDir()

    # Holds a list (by aligner) of lists (by type of statistic) of stats info
    # (that might be tuples, depending on the stat)
    # TODO: Change this to dicts by stat type.
    map_stats = []

    # make the score.results.tsv, which holds score differences and aligner/condition names.
    results_file = os.path.join(work_dir, 'score.{}.results.tsv'.format(baseline_name))
    with open(results_file, 'w') as out_results_file:
        out_results = tsv.TsvWriter(out_results_file)
        out_results.comment('diff\taligner')

        def write_tsv(comp_file, a):
            """
            Read the given comparison CSV for the given condition name, and dump
            it to the combined results file.
            """
            with open(comp_file) as comp_in:
                for line in comp_in:
                    toks = line.rstrip().split(', ')
                    out_results.line(toks[1], a)

        for name, compare_id in itertools.izip(names, compare_ids):
            compare_file = os.path.join(work_dir, '{}.compare.{}.scores'.format(name, baseline_name))
            job.fileStore.readGlobalFile(compare_id, compare_file)
            context.write_output_file(job, compare_file)
            write_tsv(compare_file, name)

            # Tabulate overall statistics
            map_stats.append([job.addChildJobFn(run_portion_worse, context, name, compare_id,
                                                cores=context.config.misc_cores, memory=context.config.misc_mem,
                                                disk=context.config.misc_disk).rv()])
            
    context.write_output_file(job, results_file)
    
    return job.addFollowOnJobFn(run_write_score_stats, context, baseline_name, names, map_stats).rv()
    
def run_write_score_stats(job, context, baseline_name, names, map_stats):
    """
    write the score comparison statistics against the baseline with the given
    name as tsv named "score.stats.<baseline name>.tsv".
    
    Returns the file ID for that file.
    
    This is different than the stats TSV format used internally, for read stats.
    """

    work_dir = job.fileStore.getLocalTempDir()
    stats_file = os.path.join(work_dir, 'score.stats.{}.tsv'.format(baseline_name))
    with open(stats_file, 'w') as stats_out_file:
        # Put each stat as a different column.
        stats_out = tsv.TsvWriter(stats_out_file)
        stats_out.comment('aligner\tcount\tworse')
        for name, stats in zip(names, map_stats):
            stats_out.line(name, stats[0][0], stats[0][1])

    return context.write_output_file(job, stats_file)
    
def run_portion_worse(job, context, name, compare_id):
    """
    Compute percentage of reads that get worse from the baseline graph.
    Return total reads and portion that got worse.
    """
    
    work_dir = job.fileStore.getLocalTempDir()

    compare_file = os.path.join(work_dir, '{}.compare.scores'.format(name))
    job.fileStore.readGlobalFile(compare_id, compare_file)
    
    total = 0
    worse = 0
    with open(compare_file) as compare_f:
        for line in compare_f:
            toks = line.split(', ')
            total += 1
            if int(toks[1]) < 0:
                worse += 1
                
    portion = float(worse) / float(total) if total > 0 else 0
    return total, portion

def run_mapeval(job, context, options, xg_file_ids, gcsa_file_ids, gbwt_file_ids, id_range_file_ids, snarl_file_ids,
                vg_file_ids, gam_file_ids, reads_gam_file_id, reads_xg_file_id, reads_bam_file_id,
                reads_fastq_file_ids,
                fasta_file_id, bwa_index_ids, bam_file_ids,
                pe_bam_file_ids, true_read_stats_file_id):
    """
    Main Toil job, and main entrypoint for use of vg_mapeval as a library.
    
    Run the analysis on the given files.
    
    TODO: Refactor to use a list of dicts/dict of listys for the indexes.
    
    Return the file IDs of result files.
    
    Returns a pair of the position comparison results and the score comparison
    results.
    
    Each result set is itself a pair, consisting of a list of per-graph file
    IDs, and an overall statistics file ID.
    
    """
    
    # This should be the only Toil job that actually uses options (in order to
    # orchestrate the right shape of workflow depending on whether we want
    # particular analyses done).
    
    # Make an indexing job
    index_job = job.addChildJobFn(run_map_eval_index,
                                  context,
                                  xg_file_ids,
                                  gcsa_file_ids,
                                  gbwt_file_ids,
                                  id_range_file_ids,
                                  snarl_file_ids,
                                  vg_file_ids,
                                  cores=context.config.misc_cores,
                                  memory=context.config.misc_mem,
                                  disk=context.config.misc_disk)

    # Extract our truth positions if no true_read_stats_file_id
    if not true_read_stats_file_id and reads_gam_file_id:
        annotate_job = index_job.addChildJobFn(annotate_gam, context, reads_xg_file_id, reads_gam_file_id)
        true_read_stats_file_id = annotate_job.addFollowOnJobFn(extract_gam_read_stats,
                                                                context, 'truth', annotate_job.rv()).rv()
    elif not true_read_stats_file_id and reads_bam_file_id:
        true_read_stats_file_id = index_job.addChildJobFn(extract_bam_read_stats,
                                                          context, 'truth', reads_bam_file_id, True).rv()

    # Extract our fastq reads so that all aligners get the exact same inputs    
    # todo: should be able to use same reads, interleaved, for both
    fq_reads_ids_bwa = reads_fastq_file_ids
    if reads_fastq_file_ids and options.bwa:
        fq_reads_ids_bwa = index_job.addChildJobFn(run_strip_fq_ext, context, reads_fastq_file_ids,
                                                   disk=context.config.alignment_disk).rv()
        
    fastq_fn = run_gam_to_fastq if reads_gam_file_id else run_bam_to_fastq
    fq_reads_ids, fq_paired_reads_ids, fq_paired_reads_for_vg_ids = (
        reads_fastq_file_ids, fq_reads_ids_bwa, reads_fastq_file_ids)
    
    # if we got two input fastqs, merge them together for single end
    if len(fq_reads_ids) == 2 and not options.paired_only:
        fq_reads_ids = [index_job.addChildJobFn(run_concat_fastqs, context, fq_reads_ids,
                                               disk=context.config.alignment_disk).rv()]
        
    if not options.paired_only and not fq_reads_ids:
        fq_reads_ids = index_job.addChildJobFn(fastq_fn, context,
                                               reads_gam_file_id if reads_gam_file_id else reads_bam_file_id,
                                               False,
                                               disk=context.config.alignment_disk).rv()
    if not options.single_only and not fq_paired_reads_ids:
        fq_paired_reads_ids  = index_job.addChildJobFn(fastq_fn, context,
                                                       reads_gam_file_id if reads_gam_file_id else reads_bam_file_id,
                                                       True,
                                                       disk=context.config.alignment_disk).rv()
        # todo: smarter annotation so we don't need to make two input sets, one with _1 _2 with vg and one without for bwa
        fq_paired_reads_for_vg_ids = index_job.addChildJobFn(fastq_fn, context,
                                                             reads_gam_file_id if reads_gam_file_id else reads_bam_file_id,
                                                             True, True,
                                                             disk=context.config.alignment_disk).rv()

    do_single_path = not options.multipath_only
    do_multi_path = options.multipath or options.multipath_only
                              
    # Then after indexing, do alignment
    alignment_job = index_job.addFollowOnJobFn(run_map_eval_align, context, index_job.rv(),
                                               options.gam_names, gam_file_ids,
                                               fq_reads_ids, fq_paired_reads_ids, fq_paired_reads_for_vg_ids,
                                               fasta_file_id, bwa_index_ids, options.bwa,
                                               not options.paired_only, not options.single_only,
                                               do_single_path, do_multi_path, 
                                               options.ignore_quals)
                                               
    # Unpack the alignment job's return values
    # TODO: we're clobbering input values...
    (gam_names, gam_file_ids, xg_ids, vg_map_times, bwa_bam_file_ids, bwa_map_times) = (
        alignment_job.rv(0), alignment_job.rv(1), alignment_job.rv(2),
        alignment_job.rv(3), alignment_job.rv(4), alignment_job.rv(5))

    # We make a root for comparison here to encapsulate its follow-on chain
    comparison_parent_job = Job()
    alignment_job.addFollowOn(comparison_parent_job)
    
    # Then do mapping evaluation comparison (the rest of the workflow)
    comparison_job = comparison_parent_job.addChildJobFn(run_map_eval_comparison, context, xg_ids,
                     gam_names, gam_file_ids, options.bam_names, bam_file_ids,
                     options.pe_bam_names, pe_bam_file_ids, bwa_bam_file_ids,
                     true_read_stats_file_id, options.mapeval_threshold, options.compare_gam_scores, reads_gam_file_id,
                     cores=context.config.misc_cores, memory=context.config.misc_mem,
                     disk=context.config.misc_disk)

    # Then do the R plotting
    position_comparison_results = comparison_job.rv(0)
    score_comparison_results= comparison_job.rv(1)
    plot_job = comparison_parent_job.addFollowOnJobFn(run_map_eval_plot, context,  position_comparison_results, score_comparison_results)

    # Dump out the running times into map_times.tsv
    comparison_parent_job.addChildJobFn(run_write_map_times, context, gam_names, vg_map_times, bwa_map_times)
                     
    return comparison_job.rv()

def run_map_eval_plot(job, context, position_comp_results, score_comp_results):
    """
    Make the PR and QQ plots with R
    """
    work_dir = job.fileStore.getLocalTempDir()

    position_stats_file_id = position_comp_results[1]
    position_stats_path = os.path.join(work_dir, 'position_stats.tsv')
    job.fileStore.readGlobalFile(position_stats_file_id, position_stats_path)

    out_name_id_pairs = []
    for rscript in ['pr', 'qq', 'roc']:

        plot_name = 'plot-{}.svg'.format(rscript)
        script_path = get_vg_script(job, context.runner, 'plot-{}.R'.format(rscript), work_dir)
        cmd = ['Rscript', os.path.basename(script_path), os.path.basename(position_stats_path),
               plot_name]
        try:
            context.runner.call(job, cmd, work_dir = work_dir)
            out_name_id_pairs.append((plot_name, context.write_output_file(job, os.path.join(work_dir, plot_name))))
        except Exception as e:
            if rscript == 'roc':
                logger.warning('plot-roc.R failed: '.format(str(e)))
            else:
                # We insist that the R scripts execute successfully (except plot-roc)
                raise e
            
    return out_name_id_pairs

def run_write_map_times(job, context, gam_names, vg_map_times, bwa_map_times):
    """
    Make a table of running times (in seconds) for mapping.  These times do not include 
    toil-vg overhead like downloading and chunking
    """

    work_dir = job.fileStore.getLocalTempDir()
    times_path = os.path.join(work_dir, 'map_times.tsv')
    with open(times_path, 'w') as times_file:
        times_file.write('aligner\tmap time (s)\n')
        for name, map_time in zip(gam_names, vg_map_times):
            if map_time is not None:
                times_file.write('{}\t{}\n'.format(name, round(map_time, 5)))
        if bwa_map_times[0] is not None:
            times_file.write('bwa-mem\t{}\n'.format(round(bwa_map_times[0], 5)))
        if bwa_map_times[1] is not None:
            times_file.write('bwa-mem-pe\t{}\n'.format(round(bwa_map_times[1], 5)))

    context.write_output_file(job, times_path)

def make_mapeval_plan(toil, options):
    """
    Import all the necessary files form options into Toil.
    
    Keep the IDs under names in an argparse namespace that functions as a "plan"
    for the workflow.
    
    """
    
    # Make a plan
    plan = argparse.Namespace()
    
    start_time = timeit.default_timer()
            
    # Upload local files to the remote IO Store
    
    # Input vg data (either pre-aligned or to align against). Can be either .vg
    # (to index and align against) or .xg/.gcsa/.gcsa.lcp (to align against) or
    # .gam/.xg (pre-alligned, just annotate)
    plan.gam_file_ids = []
    if options.gams:
        for gam in options.gams:
            plan.gam_file_ids.append(toil.importFile(gam))

    plan.vg_file_ids = []
    if options.vg_graphs:
        for graph in options.vg_graphs:
            plan.vg_file_ids.append(toil.importFile(graph))

    plan.xg_file_ids = []
    plan.gcsa_file_ids = [] # list of gcsa/lcp pairs
    plan.gbwt_file_ids = []
    plan.id_range_file_ids = []
    plan.snarl_file_ids = []
    imported_xgs = {}
    if options.index_bases:
        for ib in options.index_bases:
            imported_xgs[ib + '.xg'] = toil.importFile(ib + '.xg')
            plan.xg_file_ids.append(imported_xgs[ib + '.xg'])
            if not options.gams:
                plan.gcsa_file_ids.append(
                    (toil.importFile(ib + '.gcsa'),
                    toil.importFile(ib + '.gcsa.lcp')))
                    
                if options.use_gbwt:
                    plan.gbwt_file_ids.append(toil.importFile(ib + '.gbwt'))
                    
                if options.use_snarls:
                    plan.snarl_file_ids.append(toil.importFile(ib + '.snarls'))
                    
                # multiple gam outputs not currently supported by evaluation pipeline
                #if os.path.isfile(os.path.join(ib, '_id_ranges.tsv')):
                #    id_range_file_ids.append(
                #        toil.importFile(ib + '_id_ranges.tsv'))

    plan.reads_xg_file_id = None
    if options.gam_input_xg:
        if options.gam_input_xg in imported_xgs:
            plan.reads_xg_file_id = imported_xgs[options.gam_input_xg]
        else:
            plan.reads_xg_file_id = toil.importFile(options.gam_input_xg)
                    
    # Import input reads to be realigned
    if options.gam_input_reads:
        plan.reads_gam_file_id = toil.importFile(options.gam_input_reads)
    else:
        plan.reads_gam_file_id = None

    # Import input reads to be realigned
    if options.bam_input_reads:
        plan.reads_bam_file_id = toil.importFile(options.bam_input_reads)
    else:
        plan.reads_bam_file_id = None

    # Import input reads to be realigned        
    plan.reads_fastq_file_ids = []
    if options.fastq:
        for sample_reads in options.fastq:
            plan.reads_fastq_file_ids.append(toil.importFile(sample_reads))
                                
    # Input bwa data        
    plan.bam_file_ids = []
    if options.bams:
        for bam in options.bams:
            plan.bam_file_ids.append(toil.importFile(bam))
    plan.pe_bam_file_ids = []
    if options.pe_bams:
        for bam in options.pe_bams:
            plan.pe_bam_file_ids.append(toil.importFile(bam))
            
    plan.fasta_file_id = None
    plan.bwa_index_ids = None
    if options.fasta:
        plan.bwa_index_ids = dict()
        for suf in ['.amb', '.ann', '.bwt', '.pac', '.sa']:
            fidx = '{}{}'.format(options.fasta, suf)
            try:
                plan.bwa_index_ids[suf] = toil.importFile(fidx)
            except:
                logger.info('No bwa index found for {}, will regenerate'.format(options.fasta))
                plan.bwa_index_ids = None
                break
        if not plan.bwa_index_ids:
            plan.fasta_file_id = toil.importFile(options.fasta)
    if options.truth:
        plan.true_read_stats_file_id = toil.importFile(options.truth)
    else:
        plan.true_read_stats_file_id = None

    end_time = timeit.default_timer()
    logger.info('Imported input files into Toil in {} seconds'.format(end_time - start_time))
    
    return plan

def mapeval_main(context, options):
    """
    Run the mapeval workflow.
    """

    # Make sure the options are good
    validate_options(options)
    
    # How long did it take to run the entire pipeline, in seconds?
    run_time_pipeline = None
        
    # Mark when we start the pipeline
    start_time_pipeline = timeit.default_timer()

    t = copy.deepcopy(context)
    with context.get_toil(options.jobStore) as toil:
        if not toil.options.restart:

            # Import everything
            plan = make_mapeval_plan(toil, options)

            # Make a job to run the mapeval workflow, using all these various imported files.
            main_job = Job.wrapJobFn(run_mapeval,
                                     context, 
                                     options, 
                                     plan.xg_file_ids,
                                     plan.gcsa_file_ids,
                                     plan.gbwt_file_ids,
                                     plan.id_range_file_ids,
                                     plan.snarl_file_ids,
                                     plan.vg_file_ids, 
                                     plan.gam_file_ids,
                                     plan.reads_gam_file_id,
                                     plan.reads_xg_file_id,
                                     plan.reads_bam_file_id,
                                     plan.reads_fastq_file_ids,
                                     plan.fasta_file_id, 
                                     plan.bwa_index_ids, 
                                     plan.bam_file_ids,
                                     plan.pe_bam_file_ids, 
                                     plan.true_read_stats_file_id)
                
            # Output files all live in the out_store, but if we wanted to we could export them also/instead.

            # Init the outstore
            init_job = Job.wrapJobFn(run_write_info_to_outstore, context, sys.argv)
            init_job.addFollowOn(main_job)

            # Run the root job
            toil.start(init_job)
        else:
            toil.restart()
            
    end_time_pipeline = timeit.default_timer()
    run_time_pipeline = end_time_pipeline - start_time_pipeline
 
    print("All jobs completed successfully. Pipeline took {} seconds.".format(run_time_pipeline))
    
