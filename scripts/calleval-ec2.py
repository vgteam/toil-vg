#!/usr/bin/python2.7
"""
Do calling comparisons on a set of alignments created by mapeval
Currently needs to be run as follows: scripts/calleval-ec2.py
"""
import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("basename", help="input file prefix.  will look for [outstore]/[basename].xg/gcsa,"
                        " [basename]_primary.xg/gcsa etc. as made by construct-hs37d5-ec2.py")
    parser.add_argument("--leader", help="name of leader created with create-ec2-leader.sh")    
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--chroms", nargs='+', default=None,
                        help="chromosome(s) (default: everything in fasta)")
    parser.add_argument("--fasta", help="fasta. ideally already indexed with bwa index", required=True)
    parser.add_argument("--truth", help="vcf baseline to compare with", required=True)
    parser.add_argument("--names", nargs="*", default=["primary", "minaf_0.034"],
                        help="related graphs to try in addtion to xg_basename")
    parser.add_argument("--outname", required=True, help="subdirectory in outstore for all output (same as mapeval-ec2.py)")    
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")    
    args = args[1:]
    return parser.parse_args(args)
options = parse_args(sys.argv)

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)

s3_outstore = options.out_store.replace('aws:us-west-2:', 's3://')
me_outstore = os.path.join(options.out_store, options.outname)
s3_me_outstore = me_outstore.replace('aws:us-west-2:', 's3://')
    
log_name = 'calleval_{}.log'.format(os.path.basename(options.basename))
os_log_name = os.path.join(s3_outstore, options.outname, os.path.basename(log_name))

cmd = ['calleval', options.job_store, me_outstore,
       '--logFile', log_name,       
       '--vcfeval_fasta', options.fasta,
       '--vcfeval_baseline', options.truth,
       '--sample_name', 'SAMPLE',
       '--call_and_genotype']

if options.chroms:
    cmd += ['--chroms'] + options.chroms

# make basename absolute
absolute_basename = os.path.join(s3_outstore, options.basename)

cmd += ['--xg_paths', absolute_basename + '.xg']
for name in options.names:
    if name == 'primary':
        cmd += [os.path.join(os.path.dirname(absolute_basename), name) + '.xg']
    else:
        cmd += [absolute_basename + '_' + name + '.xg']

cmd += ['--gams']
for name in [os.path.basename(options.basename)] + options.names:
    cmd += [os.path.join(s3_me_outstore, 'aligned-{}-pe_default.gam'.format(name))]
cmd += ['--gam_names', os.path.basename(options.basename)] + options.names

cmd += ['--bams', os.path.join(s3_me_outstore, 'bwa-mem-pe.bam')]
cmd += ['--bam_names', 'bwa-mem']

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
