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
    parser.add_argument("--gbwt", action="store_true", help="make gbwt")
    parser.add_argument("--xg", action="store_true", help="make xg")
    parser.add_argument("--gcsa", action="store_true", help="make gcsa")
    parser.add_argument("--snarls", action="store_true", help="make snarls")
    parser.add_argument("--control", help="control sample")
    parser.add_argument("--primary", action="store_true", help="make primary graph")
    parser.add_argument("--minaf", type=float, help="make min allele filter graph")
    parser.add_argument("--alt_paths", action="store_true", help="force alt paths")
    parser.add_argument("--filter_ceph", action="store_true", help="filter private CEPH variants")
    
    parser.add_argument("--node", help="toil node type (default=r3.8xlarge:0.85)", default="r3.8xlarge:0.85")
    args = args[1:]        
    return parser.parse_args(args)
options = parse_args(sys.argv)

def get_vcf_path_hs37d5(chrom):
    base = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
    try:
        if int(chrom) in range(1, 23):
            return os.path.join(base, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom))
    except:
        pass
    if chrom == 'X':
        return os.path.join(base, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')
    elif chrom == 'Y':
        return os.path.join(base, 'ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz')

def get_unphased_vcf_path_hs37d5():
    return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'

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
       '--logFile', log_name]

# Note config file path is on the leader!!!!  Should fix to copy it over, but not sure how.
cmd += ['--config', options.config] if options.config else ['--whole_genome_config']

if options.minaf:
    cmd += ['--min_af', str(options.minaf)]

if options.primary:
    cmd += ['--primary']

if options.xg:
    cmd += ['--xg_index']

if options.gcsa:
    cmd += ['--gcsa_index']

if options.gbwt:
    cmd += ['--gbwt_index', '--gbwt_prune']
    
if options.gbwt or options.alt_paths:
    cmd += ['--alt_paths']

if options.snarls:
    cmd += ['--snarls_index']

if options.control:
    cmd += ['--control_sample', options.control]
    cmd += ['--filter_sample', options.control]

if options.filter_ceph:
    cmd += ['--filter_ceph']

if options.chroms:
    # restrict to specified chromosome(s)
    cmd += ['--regions'] + options.chroms
    if options.gbwt or options.alt_paths or options.filter_ceph:
        cmd += ['--vcf'] + [get_vcf_path_hs37d5(chrom) for chrom in options.chroms]
    else:
        cmd += ['--vcf'] + [get_unphased_vcf_path_hs37d5()]
else:
    # do all chromsomes as well as decoys
    cmd += ['--fasta_regions', '--vcf']
    if options.gbwt or options.alt_paths or options.filter_ceph:
        cmd += [get_vcf_path_hs37d5(chrom) for chrom in range(1, 23) + ['X', 'Y']]
    else:
        cmd += [get_unphased_vcf_path_hs37d5()]

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

