#!/usr/bin/python2.7
"""
Do mapping comparisons on a set of graphs
Currently needs to be run as follows: scripts/eval-ec2.py
"""
import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("truth_reads", help="annotated gam or bam or fastq to use as truth. if not bam, .pos file needed too")
    parser.add_argument("--leader", help="name of leader created with create-ec2-leader.sh")    
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--fasta", help="fasta. ideally already indexed with bwa index", required=True)
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")
    parser.add_argument("--outname", required=True, help="subdirectory in outstore for all output")
    parser.add_argument("--index_bases", required=True, nargs='+', help="locations of input indexes (no extension)")
    parser.add_argument("--names", required=True, nargs='+', help="names of indexes")
    args = args[1:]
    return parser.parse_args(args)
options = parse_args(sys.argv)

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)

basename = os.path.basename(options.index_bases[0])
s3_outstore = options.out_store.replace('aws:us-west-2:', 's3://')
me_outstore = os.path.join(options.out_store, options.outname)
    
log_name = 'mapeval_{}.log'.format(os.path.basename(basename))
os_log_name = os.path.join(s3_outstore, options.outname, os.path.basename(log_name))

cmd = ['mapeval', options.job_store, me_outstore,
       '--fasta', options.fasta,
       '--logFile', log_name,       
       '--bwa',
       '--multipath']

# make basename absolute
absolute_basename = os.path.join(s3_outstore, basename)

if options.truth_reads.endswith('.bam'):
    cmd += ['--bam_input_reads', options.truth_reads]
else:
    if options.truth_reads.endswith('.gam'):
        cmd += ['--gam_input_reads', options.truth_reads]
    else:
        cmd += ['--fastq', options.truth_reads]
    true_pos = options.truth_reads
    if true_pos.endswith('.fq.gz'):
        true_pos = true_pos[:-3]
    cmd += ['--truth', os.path.splitext(true_pos)[0] + '.pos']

cmd += ['--index-bases'] + options.index_bases
cmd += ['--gam-names'] + options.names

if 'ekg' in options.truth_reads:
    cmd += ['--ignore-quals']

# Note config file path is on the leader!!!!  Should fix to copy it over, but not sure how.
cmd += ['--config', options.config] if options.config else ['--whole_genome_config']

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
           '/venv/bin/aws', 's3', 'cp', log_name, os_log_name]
    print ' '.join(cmd)
    subprocess.check_call(cmd)
