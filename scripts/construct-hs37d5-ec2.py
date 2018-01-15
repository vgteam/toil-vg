#!/usr/bin/python2.7
"""
Construct set of "snp1kg" graphs and their indexes on EC2.  Wraps toil-vg construct
but fills in many input and Toil parameters to make it easier to use. 
Currently needs to be run as follows: scripts/construct-hs37d5-ec2.py
"""

import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("leader", help="name of leader created with create-ec2-leader.sh")
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("--chroms", nargs='+', default=None,
                        help="chromosome(s) (default: everything in fasta)")
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")
    parser.add_argument("--gbwt", action="store_true", help="make gbwt (only works with 1 chrom currently)")
    parser.add_argument("--control", help="control sample")
    parser.add_argument("--node", help="toil node type (default=r3.8xlarge:0.85)", default="r3.8xlarge:0.85")
    args = args[1:]        
    return parser.parse_args(args)
options = parse_args(sys.argv)

def get_vcf_path_hs37d5(chrom):
    return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom)

def get_fasta_path_hs37d5():
    return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)    
    
out_name = 'snp1kg' if not options.chroms else 'snp1kg_{}'.format('_'.join(options.chroms))
log_name = '/construct_{}.log'.format(out_name)
os_log_name = os.path.join(options.out_store[options.out_store.rfind(':')+1:], os.path.basename(log_name))

cmd = ['construct', options.job_store, options.out_store,
       '--fasta', get_fasta_path_hs37d5(),
       '--out_name', out_name,
       '--logFile', log_name,
       '--primary',
       '--alt_paths',
       '--min_af', '0.034',
       '--xg_index', '--gcsa_index']

# Note config file path is on the leader!!!!  Should fix to copy it over, but not sure how.
cmd += ['--config', options.config] if options.config else ['--whole_genome_config']

if options.gbwt:
    cmd += ['--gbwt_index']

if options.control:
    cmd += ['--control_sample', options.control]
    cmd += ['--filter_sample', options.control]

if options.chroms:
    # restrict to specified chromosome(s)
    cmd += ['--regions'] + options.chroms
    cmd += ['--vcf'] + [get_vcf_path_hs37d5(chrom) for chrom in options.chroms]
else:
    # do all chromsomes as well as decoys
    cmd += ['--fasta_regions']
    cmd += ['--vcf'] + [get_vcf_path_hs37d5(chrom) for chrom in range(1, 23) + ['X', 'Y']]

if options.restart:
    cmd += ['--restart']
else:
    subprocess.check_call(['toil', 'clean', options.job_store])
    
print ' '.join(cmd)
subprocess.check_call(['scripts/ec2-run.sh', '-n', options.node, options.leader, ' '.join(cmd)])

#copy the log to the out store
cmd = ['toil', 'ssh-cluster',  '--insecure', '--zone=us-west-2a', options.leader,
       '/venv/bin/aws', 's3', 'cp', log_name, 's3://{}'.format(os_log_name)]
print ' '.join(cmd)
subprocess.check_call(cmd)

