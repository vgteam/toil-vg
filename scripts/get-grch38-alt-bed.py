#!/usr/bin/python2.7
"""
Make a BED file mapping alt contigs to the positions on the reference.  I'm sure such
a file already exists somewhere, but I've yet to find it.
"""

import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fai", default='ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai',
                        help='fai index used for sequence lengths')
    parser.add_argument('--bwakit', default='https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.15_x64-linux.tar.bz2/download', help='bwa-kit tarball to fish mapping out of')
    args = args[1:]        
    return parser.parse_args(args)

options = parse_args(sys.argv)

# download
subprocess.check_call(['wget', '-nc', options.fai])
subprocess.check_call(['wget', '-nc', options.bwakit])
subprocess.check_call(['tar', 'xjf', os.path.basename(options.bwakit)])
alts_path = 'bwa.kit/resource-GRCh38/hs38DH.fa.alt'

# get our lengths from the fai (note these could also be parsed from the cigar strings in alts_path)
with open(os.path.basename(options.fai)) as fai_file:
    len_map = {}
    for line in fai_file:
        toks = line.strip().split('\t')
        if len(toks) >=2 and not toks[0].startswith('#'):
            len_map[toks[0]] = int(toks[1])

# spit out our bed
with open(alts_path) as alts_file:
    for line in alts_file:
        toks = line.strip().split('\t')
        if len(toks) >= 3 and not toks[0].startswith('@') and int(toks[3]) != 0:
            sys.stdout.write('{}\t{}\t{}\t{}\n'.format(toks[2], int(toks[3]) - 1,
                                                       int(toks[3]) + len_map[toks[0]], toks[0]))

sys.stderr.write('\n*****\nYou may want to run\n\n rm -rf {} {} {}\n\nto clean up\n*****\n'.format(
    'bwa-kit', os.path.basename(options.fai), os.path.basename(options.bwakit)))
            
