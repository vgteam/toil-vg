#!/usr/bin/python2.7
"""
Construct a baseline graph to simulate from for testing. Wraps toil-vg construct
but fills in many input and Toil parameters to make it easier to use. 
Currently needs to be run as follows: scripts/construct-hs37d5-baseline-ec2.py
"""

import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("--leader", help="name of leader created with create-ec2-leader.sh")    
    parser.add_argument("--chroms", nargs='+', default=[str(x) for x in range(1,23)],
                        help="chromosome(s) (default: 1-22)")
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")
    parser.add_argument("--control", help="control sample", required=True)
    args = args[1:]        
    return parser.parse_args(args)
options = parse_args(sys.argv)

def get_vcf_paths_hs37d5(chroms, sample):
    if sample == 'HG002':
        return ['ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz']
    else:
        return ['ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom) for chrom in chroms]

def get_fasta_path_hs37d5():
    return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)    
    
out_name = 'baseline'
if options.chroms != [str(x) for x in range(1,23)]:
    out_name += '_' + '_'.join(options.chroms)
log_name = 'construct_{}.log'.format(out_name)
os_log_name = os.path.join(options.out_store[options.out_store.rfind(':')+1:], os.path.basename(log_name))

cmd = ['construct', options.job_store, options.out_store,
       '--fasta', get_fasta_path_hs37d5(),
       '--out_name', out_name,
       '--logFile', log_name,
       '--haplo_sample', options.control,
       '--control_sample', options.control,
       '--alt_paths',
       '--xg_index']

# Note config file path is on the leader!!!!  Should fix to copy it over, but not sure how.
cmd += ['--config', options.config] if options.config else ['--whole_genome_config']
cmd += ['--control_sample', options.control]
cmd += ['--regions'] + options.chroms
cmd += ['--vcf'] + get_vcf_paths_hs37d5(options.chroms, options.control)

if options.restart:
    cmd += ['--restart']
else:
    subprocess.check_call(['toil', 'clean', options.job_store])

ec2_cmd = ['scripts/ec2-run.sh']
if options.leader:
    ec2_cmd += ['-l', options.leader]
ec2_cmd += [' '.join(cmd)]

print ' '.join(ec2_cmd)
subprocess.check_call(ec2_cmd)

#copy the log to the out store
if options.leader:
    cmd = ['toil', 'ssh-cluster',  '--insecure', '--zone=us-west-2a', options.leader,
           '/venv/bin/aws', 's3', 'cp', log_name, 's3://{}'.format(os_log_name)]
    print ' '.join(cmd)
    subprocess.check_call(cmd)

