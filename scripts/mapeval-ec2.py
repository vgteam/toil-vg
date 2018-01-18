#!/usr/bin/python2.7
"""
Do mapping comparisons on a set of graphs
Currently needs to be run as follows: scripts/eval-ec2.py
"""
import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("leader", help="name of leader created with create-ec2-leader.sh")
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("basename", help="input file prefix.  will look for [outstore]/[basename].xg/gcsa,"
                        " [basename]_primary.xg/gcsa etc. as made by construct-hs37d5-ec2.py")
    parser.add_argument("truth_reads", help="annotated gam or bam to use as truth. if gam, .pos file needed too")
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--fasta", help="fasta. ideally already indexed with bwa index", required=True)
    parser.add_argument("--names", nargs="*", default=["primary", "minaf_0.034"],
                        help="related graphs to try in addtion to basename"),
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")
    parser.add_argument("--outname", required=True, help="subdirectory in outstore for all output")
    args = args[1:]
    return parser.parse_args(args)
options = parse_args(sys.argv)

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)

s3_outstore = options.out_store.replace('aws:us-west-2:', 's3://')
me_outstore = os.path.join(options.out_store, options.outname)
    
log_name = '/mapeval_{}.log'.format(os.path.basename(options.basename))
os_log_name = os.path.join(s3_outstore, options.outname, os.path.basename(log_name))

cmd = ['mapeval', options.job_store, me_outstore,
       '--fasta', options.fasta,
       '--logFile', log_name,       
       '--bwa',
       '--multipath']

# make basename absolute
absolute_basename = os.path.join(s3_outstore, options.basename)
absolute_truth = os.path.join(s3_outstore, options.truth_reads)

if options.truth_reads.endswith('.bam'):
    cmd += ['--bam_input_reads', absolute_truth]
else:
    cmd += ['--gam_input_reads', absolute_truth]
    cmd += ['--truth', os.path.splitext(absolute_truth)[0] + '.pos']

cmd += ['--index-bases', absolute_basename]
for name in options.names:
    if name == 'primary':
        cmd += [os.path.join(os.path.dirname(absolute_basename), name)]
    else:
        cmd += [absolute_basename + '_' + name]

cmd += ['--gam-names', os.path.basename(options.basename)] + options.names

if 'ekg' in options.truth_reads:
    cmd += ['--ignore-quals']

# Note config file path is on the leader!!!!  Should fix to copy it over, but not sure how.
cmd += ['--config', options.config] if options.config else ['--whole_genome_config']

if options.restart:
    cmd += ['--restart']
else:
    subprocess.check_call(['toil', 'clean', options.job_store])
    
print ' '.join(cmd)
subprocess.check_call(['scripts/ec2-run.sh', options.leader, ' '.join(cmd)])

#copy the log to the out store
cmd = ['toil', 'ssh-cluster',  '--insecure', '--zone=us-west-2a', options.leader,
       '/venv/bin/aws', 's3', 'cp', log_name, os_log_name]
print ' '.join(cmd)
subprocess.check_call(cmd)
