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

        # Optional: Use the following reference paths.  Some possibilities for whole human genome below:
        # path-name: ['1', '2', '3', '4', '5', '6', '7', '8' '9' '10', '11', '12', '13', '14', '15', 16', '17', '18', '19', '20', '21', '22', 'X', 'Y']
        # path-name: ['1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8' '9' '10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']

        
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
        
        # Core arguments for vg mapping
        # Note -i/--interleaved will be ignored. use the --interleaved option 
        # on the toil-vg command line instead
        vg-map-args: ['-M2', '-W', '500', '-u', '0', '-U', '-O', '-S', '50', '-a']
        
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


def apply_config_file_args(args):
    """
    Merge args from the config file and the parser, giving priority to the parser.
    """

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
