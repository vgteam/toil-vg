#!/usr/bin/env python2.7
"""
vg_config.py: Default configuration values all here (and only here), as well as logic
for reading and generating config files.

"""
from __future__ import print_function
import argparse, sys, os, os.path, errno, random, subprocess, shutil, itertools, glob, tarfile
import doctest, re, json, collections, time, timeit
import logging, logging.handlers, SocketServer, struct, socket, threading
import string
import urlparse
import getpass
import pdb
import textwrap
import yaml
from toil_lib import require

# TODO: configure RPATH-equivalent on OS X for finding libraries without environment variables at runtime
def generate_config():
    return textwrap.dedent("""
        # Toil VG Pipeline configuration file (created by toil-vg generate-config)
        # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
        # Edit the values in the configuration file and then rerun the pipeline: "toil-vg run"
        # 
        # URLs can take the form: "/", "s3://"
        # Local inputs follow the URL convention: "/full/path/to/input.txt"
        # S3 URLs follow the convention: "s3://bucket/directory/file.txt"
        #
        # Comments (beginning with #) do not need to be removed. 
        # Command-line options take priority over parameters in this file.  
        ######################################################################################################################

        ###########################################
        ### Toil resource tuning                ###
        
        # These parameters must be adjusted based on data and cluster size
        # when running on anything other than single-machine mode

        # TODO: Reduce number of parameters here.  Seems fine grained, especially for disk/mem
        # option to spin off config files for small/medium/large datasets?   

        # The following parameters assign resources to small helper jobs that typically don't do 
        # do any computing outside of toil overhead.  Generally do not need to be changed. 
        misc-cores: 1
        misc-mem: '1G'
        misc-disk: '1G'

        # Resources allotted for xg and gcaa indexing. 
        index-cores: 1
        index-mem: '4G'
        index-disk: '2G'

        # Resources for fastq splitting and gam merging
        # fastq spilt can use up to 2 cores, gam merging single threaded
        fq-split-cores: 1
        fq-split-mem: '4G'
        fq-split-disk: '2G'

        # Number of threads to use for Rocksdb GAM indexing
        # Generally, this should be kept low as speedup drops off radically 
        # after a few threads.
        gam-index-cores: 1

        # Resources for *each* vg map job
        # the number of vg map jobs is controlled by reads-per-chunk (below)
        alignment-cores: 1
        alignment-mem: '4G'
        alignment-disk: '2G'

        # Resources for chunking up a graph/gam for calling (and merging)
        # typically take xg for whoe grpah, and gam for a chromosome,
        # and split up into chunks of call-chunk-size (below)
        call-chunk-cores: 1
        call-chunk-mem: '4G'
        call-chunk-disk: '2G'

        # Resources for calling each chunk (currently includes pileup/call/genotype)
        calling-cores: 1
        calling-mem: '4G'
        calling-disk: '2G'

        ###########################################
        ### Arguments Shared Between Components ###
        # Use output store instead of toil for all intermediate files (use only for debugging)
        force-outstore: False

        # Use the following reference paths.  Some possibilities for whole human genome below:
        # path-name: ['1', '2', '3', '4', '5', '6', '7', '8' '9' '10', '11', '12', '13', '14', '15', 16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
        # path-name: ['1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8' '9' '10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
        
        #############################
        ### Docker Tool Arguments ###
        # Do not use docker for any commands
        no-docker: False
        
        ## Docker Tool List ##
        ##   Each tool is specified as a list where the first element is the docker image URL,
        ##   and the second element indicates if the docker image has an entrypoint or not
        ##   If left blank or commented, then the tool will be run directly from the command line instead
        ##   of through docker. no-docker (above) overrides all these options. 
        # Dockerfile to use for vg

        vg-docker: ['quay.io/glennhickey/vg:v1.4.0-1980-g38453bf', True]

        # Dockerfile to use for bcftools
        bcftools-docker: ['quay.io/cmarkello/bcftools', False]

        # Dockerfile to use for tabix
        tabix-docker: ['quay.io/cmarkello/htslib:latest', False]

        # Dockerfile to use for jq
        jq-docker: ['devorbitus/ubuntu-bash-jq-curl', False]
        
        # Dockerfile to use for rtg
        rtg-docker: ['realtimegenomics/rtg-tools:3.7.1', True]
        
        ##########################
        ### vg_index Arguments ###
        # Maximum edges to cross in index
        edge-max: 5

        # Size of kmers to use in indexing and mapping
        kmer-size: 10

        # Use the pruned graph in the index
        include-pruned: False

        ########################
        ### vg_map Arguments ###

        # Number of reads per chunk to use when splitting up fastq.  
        # Each chunk will correspond to a vg map job
        reads-per-chunk: 10000000
        
        # Context expansion used for gam chunking
        chunk_context: 50
    
        # Treat input fastq as paired-end interleaved
        interleaved: False
        
        # Core arguments for vg mapping (do not include file names or -t/--threads)
        # Note -i/--interleaved will be ignored. use the --interleaved option 
        # on the toil-vg command line instead
        map-opts: ['-M2', '-W', '500', '-u', '0', '-U', '-O', '-S', '50', '-a']
        
        # Type of vg index to use for mapping (either 'gcsa-kmer' or 'gcsa-mem')
        index-mode: gcsa-mem

        #########################
        ### vg_call Arguments ###
        # Overlap option that is passed into make_chunks and call_chunk
        overlap: 2000
        
        # Chunk size
        call-chunk-size: 10000000

        # Options to pass to chunk_gam. (do not include file names or -t/--threads)
        filter-opts: ['-r', '0.9', '-fu', '-s', '1000', '-o', '0', '-q', '15']

        # Options to pass to vg pileup. (do not include file names or -t/--threads)
        pileup-opts: ['-q', '10', '-a']

        # Options to pass to vg call. (do not include file names or -t/--threads)
        call-opts: ['']
        
        # Options to pass to vg genotype. (do not include file names or -t/--threads)
        genotype-opts: ['']

        # Use vg genotype instead of vg call
        genotype: False
        
    """)


def apply_config_file_args(args):
    """
    Merge args from the config file and the parser, giving priority to the parser.
    """

    # turn --*_opts from strings to lists to be consistent with config file
    for x_opts in ['map_opts', 'call_opts', 'filter_opts', 'genotype_opts']:
        if x_opts in args.__dict__.keys() and type(args.__dict__[x_opts]) is str:
            args.__dict__[x_opts] = args.__dict__[x_opts].split(' ')
            # get rid of any -t or --threads while we're at it
            for t in ['-t', '--threads']:
                if t in args.__dict__[x_opts]:
                    pos = args.__dict__[x_opts].index(t)
                    del args.__dict__[x_opts][pos:pos+2]

    # If no config file given, we generate a default one
    if args.config is None:
        config = generate_config()
    else:
        require(os.path.exists(args.config), 'Config, {}, not found. Please run '
            '"toil-vg generate-config > {}" to create.'.format(args.config, args.config))    
        with open(args.config) as conf:
            config = conf.read()
                
    # Parse config
    parsed_config = {x.replace('-', '_'): y for x, y in yaml.load(config).iteritems()}
    options = argparse.Namespace(**parsed_config)

    # Add in options from the program arguments to the arguments in the config file
    #   program arguments that are also present in the config file will overwrite the
    #   arguments in the config file
    for args_key in args.__dict__:
        # Add in missing program arguments to config option list and
        # overwrite config options with corresponding options that are not None in program arguments
        if (args.__dict__[args_key]) or (args_key not in  options.__dict__.keys()):
            options.__dict__[args_key] = args.__dict__[args_key]
            
    return options

def config_main(args):
    """ config just prints out a file """
    
    sys.stdout.write(generate_config())
