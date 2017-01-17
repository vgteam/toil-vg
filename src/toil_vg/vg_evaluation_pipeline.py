#!/usr/bin/env python2.7
"""
vg_evaluation_pipeline.py: Run the mapping and variant calling evaluation on all the servers in
parallel using the vg framework.

old docker vg tool=1.4.0--4cbd3aa6d2c0449730975517fc542775f74910f3
new docker vg tool=latest

chr_length_list obtained from Mike Lin's vg dnanexus pipeline configuration
    chr_label_list = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"]
    chr_length_list = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import getpass
import textwrap
import yaml

from math import ceil
from subprocess import Popen, PIPE
from Bio import SeqIO

from toil.common import Toil
from toil.job import Job
from toil_lib import require
from toil_lib.toillib import *
from toil_vg.vg_common import *
from toil_vg.vg_call import *
from toil_vg.vg_index import *
from toil_vg.vg_map import *
from toil_vg.vg_vcfeval import *

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
    parser = argparse.ArgumentParser(prog='vg_evaluation_pipeline', description=main.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    subparsers = parser.add_subparsers(dest='command')
    
    # Config subparser
    parser_config = subparsers.add_parser('generate-config',
                                          help='Generates an editable default config file')
    parser_config.add_argument('--config', default='config-toil-vg.tsv', type=str,
                               help='Path to the config file to generate. '
                               '\nDefault value: "%(default)s"')
    
    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the Toil VG DNA-seq pipeline')
    pipeline_subparser(parser_run)

    # Index subparser
    parser_index = subparsers.add_parser('index', help='Runs only the vg indexing')
    index_subparser(parser_index)

    # Map subparser
    parser_map = subparsers.add_parser('map', help='Runs only the vg mapping')
    map_subparser(parser_map)

    # Call subparser
    parser_call = subparsers.add_parser('call', help='Runs only the vg calling')
    call_subparser(parser_call)

    # vcfeval subparser
    parser_vcfeval = subparsers.add_parser('vcfeval', help='Compare two VCFs wit rtg vcfeval')
    vcfeval_subparser(parser_vcfeval)

    return parser.parse_args()


def pipeline_subparser(parser_run):
    
    # Add the Toil options so the job store is the first argument
    Job.Runner.addToilOptions(parser_run)
    
    # General options
    parser_run.add_argument("vg_graph", type=str,
        help="Input vg graph file path")
    parser_run.add_argument("sample_reads", type=str,
        help="Path to sample reads in fastq format")
    parser_run.add_argument("sample_name", type=str,
        help="sample name (ex NA12878)")
    parser_run.add_argument("out_store",
        help="output IOStore to create and fill with files that will be downloaded to the local machine where this toil script was run")
    parser_run.add_argument("--xg_index", type=str,
        help="Path to xg index (to use instead of generating new one)")    
    parser_run.add_argument("--gcsa_index", type=str,
        help="Path to GCSA index (to use instead of generating new one)")
    parser_run.add_argument("--path_name", nargs='+', type=str,
        help="Name of reference path in the graph (eg. ref or 17)")
    parser_run.add_argument("--path_size", nargs='+', type=int,
        help="Size of the reference path in the graph")    
    parser_run.add_argument("--restat", default=False, action="store_true",
        help="recompute and overwrite existing stats files")

    # Add common options shared with everybody
    add_common_vg_parse_args(parser_run)
    
    # Add common indexing options shared with vg_index
    index_parse_args(parser_run)

    # add common mapping options shared with vg_map
    map_parse_args(parser_run)
    
    # Add common calling options shared with vg_call
    chunked_call_parse_args(parser_run)

    # Add common docker options
    add_docker_tool_parse_args(parser_run)



# Below are the top level jobs of the vg_evaluation pipeline.  They
# form a chain of "follow-on" jobs as follows:
#
# run_pipeline_upload --> run_pipeline_index --> run_pipeline_map --> run_pipeline_call
# --> run_pipeline_merge_vcf
#
# Each part of this chain spawns a child job that can also be run independently
# using one of vg_index.main, vg_map.main etc.
#
# Data is communicated across the chain via the output store (at least for now). 


def run_pipeline_index(job, options, inputGraphFileID, inputReadsFileID, inputXGFileID,
                       inputGCSAFileID, inputLCPFileID):
    """
    All indexing.  result is a tarball in thie output store.
    """

    if inputXGFileID is None:
        xg_file_id = job.addChildJobFn(run_xg_indexing, options, inputGraphFileID,
                                       cores=options.index_cores, memory="4G", disk="2G").rv()
    else:
        xg_file_id = inputXGFileID
        
    if inputGCSAFileID is None:
        gcsa_and_lcp_ids = job.addChildJobFn(run_gcsa_indexing, options, inputGraphFileID,
                                       cores=options.index_cores, memory="4G", disk="2G").rv()
    else:
        assert inputLCPFileID is not None
        gcsa_and_lcp_ids = inputGCSAFileID, inputLCPFileID

    return job.addFollowOnJobFn(run_pipeline_map, options, inputGraphFileID, xg_file_id, gcsa_and_lcp_ids,
                                inputReadsFileID, cores=2, memory="4G", disk="2G").rv()

def run_pipeline_map(job, options, graph_file_id, xg_file_id, gcsa_and_lcp_ids, reads_file_id):
    """ All mapping, including fastq splitting and gam merging"""

    chr_gam_ids = job.addChildJobFn(run_split_fastq, options, graph_file_id, xg_file_id, gcsa_and_lcp_ids,
                                     reads_file_id, cores=3, memory="4G", disk="2G").rv()

    return job.addFollowOnJobFn(run_pipeline_call, options, graph_file_id, xg_file_id,
                                chr_gam_ids, cores=2, memory="4G", disk="2G").rv()

def run_pipeline_call(job, options, graph_file_id, xg_file_id, chr_gam_ids):
    """ Run variant calling on the chromosomes in parallel """

    return_value = []
    vcf_tbi_file_id_pair_list = [] 
    if options.path_name and options.path_size:
        assert len(chr_gam_ids) == len(options.path_name)
        
        #Run variant calling
        for i in range(len(chr_gam_ids)):
            alignment_file_id = chr_gam_ids[i]
            chr_label = options.path_name[i]
            chr_length = options.path_size[i]
            vcf_tbi_file_id_pair = job.addChildJobFn(run_calling, options, xg_file_id,
                                                     alignment_file_id, chr_label, chr_length,
                                                     cores=options.calling_cores, memory="4G", disk="2G").rv()
            vcf_tbi_file_id_pair_list.append(vcf_tbi_file_id_pair)
    else:
        raise RuntimeError("Invalid or non-existant path_name(s) and/or path_size(s): {}, {}".format(path_name, path_size))

    return job.addFollowOnJobFn(run_pipeline_merge_vcf, options, vcf_tbi_file_id_pair_list,
                                cores=2, memory="4G", disk="2G").rv()

def run_pipeline_merge_vcf(job, options, vcf_tbi_file_id_pair_list):

    RealTimeLogger.get().info("Completed gam merging and gam path variant calling.")
    RealTimeLogger.get().info("Starting vcf merging vcf files.")
    # Set up the IO stores each time, since we can't unpickle them on Azure for
    # some reason.
    out_store = IOStore.get(options.out_store)

    # Define work directory for docker calls
    work_dir = job.fileStore.getLocalTempDir()
    
    vcf_merging_file_key_list = [] 
    for i, vcf_tbi_file_id_pair in enumerate(vcf_tbi_file_id_pair_list):
        vcf_file = "{}/{}.gz".format(work_dir, 'vcf_chunk_{}.vcf.gz'.format(i))
        vcf_file_idx = "{}.tbi".format(vcf_file)
        read_from_store(job, options, vcf_tbi_file_id_pair[0], vcf_file)
        read_from_store(job, options, vcf_tbi_file_id_pair[1], vcf_file_idx)
        vcf_merging_file_key_list.append(os.path.basename(vcf_file))

    vcf_merged_file_key = "" 
    if len(vcf_merging_file_key_list) > 1:
        # merge vcf files
        vcf_merged_file_key = "{}.vcf.gz".format(options.sample_name)
        command=['bcftools', 'concat', '-O', 'z', '-o', os.path.basename(vcf_merged_file_key), ' '.join(vcf_merging_file_key_list)]
        options.drunner.call(command, work_dir=work_dir)
        command=['bcftools', 'tabix', '-f', '-p', 'vcf', os.path.basename(vcf_merged_file_key)]
        options.drunner.call(command, work_dir=work_dir)
    else:
        vcf_merged_file_key = vcf_merging_file_key_list[0]

    # save variant calling results to the output store
    out_store_key = "{}.vcf.gz".format(options.sample_name)
    vcf_file = os.path.join(work_dir, vcf_merged_file_key)
    vcf_file_idx = vcf_file + ".tbi"
    write_to_store(job, options, vcf_file, use_out_store = True, out_store_key = out_store_key)
    write_to_store(job, options, vcf_file_idx, use_out_store = True, out_store_key = out_store_key + ".tbi") 

def generate_config():
    return textwrap.dedent("""
        # Toil VG Pipeline configuration file
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in the configuration file and then rerun the pipeline: "toil-vg run"
        # 
        # URLs can take the form: "/", "s3://"
        # Local inputs follow the URL convention: "/full/path/to/input.txt"
        # S3 URLs follow the convention: "s3://bucket/directory/file.txt"
        #
        # Comments (beginning with #) do not need to be removed. Optional parameters left blank are treated as false or blank.
        ######################################################################################################################

        ###########################################
        ### Arguments Shared Between Components ###
        # Optional: Use output store instead of toil for all intermediate files (use only for debugging)
        force-outstore: False
        
        #############################
        ### Docker Tool Arguments ###
        # Optional: Do not use docker for any commands
        no-docker: False
        
        ## Docker Tool List ##
        ##   Each tool is specified as a list where the first element is the docker image URL,
        ##   and the second element indicates if the docker image has an entrypoint or not
        ##   If left blank, then the tool will be run directly from the command line instead
        ##   of through docker
        # Optional: Dockerfile to use for vg

        vg-docker: ['quay.io/glennhickey/vg:v1.4.0-1980-g38453bf', True]

        # Optional: Dockerfile to use for bcftools
        bcftools-docker: ['quay.io/cmarkello/bcftools', False]

        # Optional: Dockerfile to use for tabix
        tabix-docker: ['quay.io/cmarkello/htslib:latest', False]

        # Optional: Dockerfile to use for jq
        jq-docker: ['devorbitus/ubuntu-bash-jq-curl', False]
        
        # Optional: Dockerfile to use for rtg
        rtg-docker: ['realtimegenomics/rtg-tools:3.7.1', True]
        
        ##########################
        ### vg_index Arguments ###
        # Optional: Maximum edges to cross in index
        edge-max: 5

        # Optional: Size of kmers to use in indexing and mapping
        kmer-size: 10

        # Optional: Don't re-use existing indexed graphs
        reindex: False

        # Optional: Use the pruned graph in the index
        include-pruned: False

        # Optional: Use the primary path in the index
        include-primary: True

        # Optional: Number of threads during the indexing step
        index-cores: 4

        # Optional: Toil job memory allocation for indexing
        index-mem: '4G'

        # Optional: Toil job disk allocation for indexing
        index-disk: '2G'

        ########################
        ### vg_map Arguments ###
        # Optional: Number of chunks to split the input fastq file records
        num-fastq-chunks: 3
        
        # Optional: Number of threads during the alignment step
        alignment-cores: 3

        # Optional: Number of threads to use for Rocksdb GAM indexing
        # Generally, this should be kept low as speedup drops off radically 
        # after a few threads.
        gam-index-cores: 3

        # Optional: Context expansion used for gam chunking
        chunk_context: 20
        
        # Optional: Core arguments for vg mapping
        vg-map-args: ['-i', '-M2', '-W', '500', '-u', '0', '-U', '-O', '-S', '50', '-a']
        
        # Optional: Toil job memory allocation for mapping
        alignment-mem: '4G'
        
        # Optional: Toil job disk allocation for mapping
        alignment-disk: '2G'

        # Optional: Type of vg index to use for mapping (either 'gcsa-kmer' or 'gcsa-mem')
        index-mode: gcsa-mem


        #########################
        ### vg_call Arguments ###
        # Optional: Overlap option that is passed into make_chunks and call_chunk
        overlap: 2000
        
        # Optional: Chunk size
        call-chunk-size: 10000000

        # Optional: Chromosomal position offset (eg. 43044293)
        offset: None

        # Optional: Options to pass to chunk_gam.
        filter-opts: ['-r', '0.9', '-fu', '-s', '1000', '-o', '0', '-q', '15']

        # Optional: Options to pass to vg pileup.
        pileup-opts: ['-q', '10', '-a']

        # Optional: Options to pass to vg call.
        call-opts: ['']
        
        # Optional: Options to pass to vg genotype.
        genotype-opts: ['']

        # Optional: Use vg genotype instead of vg call
        genotype: False

        # Optional: Number of threads during the variant calling step
        calling-cores: 4
        
        # Optional: Toil job memory allocation for variant calling
        calling-mem: '4G'
        
        # Optional: Toil job disk allocation for variant calling
        calling-disk: '2G'
        
        # Optional: always overwrite existing files
        overwrite: False
    """)

def generate_file(file_path, generate_func):
    """
    Checks file existance, generates file, and provides message
    :param str file_path: File location to generate file
    :param function generate_func: Function used to generate file
    """
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))


def main():
    """
    Computational Genomics Lab, Genomics Institute, UC Santa Cruz
    Toil vg DNA-seq pipeline
    
    DNA-seq fastqs are split, aligned to an indexed vg reference graph, and variant-called using
    the vg toolset.

    General usage:
    Type "toil-vg [toil options] [toil-vg flag arguments] [jobStore] [path to vg file] [path to fastq file]
                  [sample name] [local output directory] [remote input fileStore] [remote output fileStore]'

    Please read the README.md located in the source directory for more documentation

    Structure of vg DNA-Seq Pipeline (per sample)
                  
                  
                       > 3 ---->    > 5 ---->
                      / ..      |  / ..     |
                     /  ..      v /  ..     v
        0 --> 1 --> 2 --> 3 --> 4 --> 5 --> 6 --> 7
                     \  ..      ^ \  ..     ^
                      \ ..      |  \ ..     | 
                       > 3 ---->    > 5 --->

    0 = Upload local input files to remote input fileStore
    1 = Index input vg reference graph
    2 = Shard Reads
    3 = Align reads to reference graph
    4 = Merge read alignments
    5 = Run vg variant calls on chromosome-chunked .gam alignment
    6 = Merge variant call files
    7 = Download output files from remote output fileStore to local output directory
    ================================================================================
    Dependencies
    toil:           pip install toil (version >= 3.5.0a1.dev241)
    toil-lib:       git clone --recursive https://github.com/cmarkello/toil-lib.git ${PWD}/toil-lib/
                    pip install ${PWD}/toil-lib/
    biopython:      pip install biopython (version >= 1.67)
    boto:           pip install boto (OPTIONAL for runs on AWS machines)
    """

    args = parse_args()
    
    # Write out our config file that's necessary for all other subcommands
    if args.command == 'generate-config':
        yaml_config_out = generate_config()
        with open(args.config, 'w') as tgt_config_file:
            tgt_config_file.write(yaml_config_out)
        return

    # Else we merge in our config file args with the command line args for
    # whatever subparser we're in
    require(os.path.exists(args.config), '{} not found Please run '
            '"toil-vg generate-config"'.format(args.config))
    # Parse config
    parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}
    options = argparse.Namespace(**parsed_config)

    # Add in options from the program arguments to the arguments in the config file
    #   program arguments that are also present in the config file will overwrite the
    #   arguments in the config file
    for args_key in args.__dict__:
        # Add in missing program arguments to config option list and
        # overwrite config options with corresponding options that are not None in program arguments
        if (args.__dict__[args_key]) or (args_key not in  options.__dict__.keys()):
            options.__dict__[args_key] = args.__dict__[args_key]
            ## Options sanity checks
            # path_name and path_size lists must be equal in length

    print('Program Argument List: ', args)
    print('Final Options List: ')
    for key in sorted(options.__dict__.keys()):
        print("option {}: {}".format(key, options.__dict__[key]))
        # If no arguments provided, print the full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # Relative outstore paths can end up who-knows-where.  Make absolute.
    if options.out_store[0] == '.':
        options.out_store = os.path.abspath(options.out_store)

    if args.command == 'run':
        pipeline_main(options)
    elif args.command == 'index':
        index_main(options)
    elif args.command == 'map':
        map_main(options)
    elif args.command == 'call':
        call_main(options)
    elif args.command == 'vcfeval':
        vcfeval_main(options)
        
    
def pipeline_main(options):
    
    RealTimeLogger.start_master()

    # make the docker runner
    options.drunner = DockerRunner(
        docker_tool_map = get_docker_tool_map(options))

    # Some file io is dependent on knowing if we're in the pipeline
    # or standalone. Hack this in here for now
    options.tool = 'pipeline'

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
            inputReadsFileID = import_to_store(toil, options, options.sample_reads)
            if options.gcsa_index:
                inputGCSAFileID = import_to_store(toil, options, options.gcsa_index)
                inputLCPFileID = import_to_store(toil, options, options.gcsa_index + ".lcp")
            else:
                inputGCSAFileID = None
                inputLCPFileID = None
            if options.xg_index:
                inputXGFileID = import_to_store(toil, options, options.xg_index)
            else:
                inputXGFileID = None

            # Make a root job
            root_job = Job.wrapJobFn(run_pipeline_index, options, inputGraphFileID,
                                     inputReadsFileID, inputXGFileID, inputGCSAFileID,
                                     inputLCPFileID,
                                     cores=2, memory="5G", disk="2G")

            # Run the job and store
            toil.start(root_job)
        else:
            toil.restart()

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

