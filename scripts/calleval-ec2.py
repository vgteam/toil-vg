#!/usr/bin/python2.7
"""
Do calling comparisons on a set of alignments created by mapeval
Currently needs to be run as follows: scripts/calleval-ec2.py
"""
import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("leader", help="name of leader created with create-ec2-leader.sh")
    parser.add_argument("job_store")
    parser.add_argument("out_store", help="must be same as used in mapeval-ec2.py")
    parser.add_argument("xg_basename", help="input file prefix.  will look for [xg_basename].xg,"
                        " [xg_basename]_primary.xg etc. as made by construct-hs37d5-ec2.py/mapeval-ec2.py")
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--chroms", nargs='+', default=None,
                        help="chromosome(s) (default: everything in fasta)")
    parser.add_argument("--fasta", help="fasta. ideally already indexed with bwa index", required=True)
    parser.add_argument("--truth", help="vcf baseline to compare with", required=True)
    parser.add_argument("--names", nargs="*", default=["primary", "minaf_0.034"],
                        help="related graphs to try in addtion to xg_basename")    
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")    
    args = args[1:]
    return parser.parse_args(args)
options = parse_args(sys.argv)

if not options.job_store.startswith('aws:'):
    options.jo_bstore = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)
    
log_name = '/calleval_{}.log'.format(os.path.basename(options.xg_basename))
os_log_name = os.path.join(options.out_store[options.out_store.rfind(':')+1:], os.path.basename(log_name))

cmd = ['calleval', options.job_store, options.out_store,
       '--vcfeval_fasta', options.fasta,
       '--vcfeval_baseline', options.truth,
       '--sample_name', 'SAMPLE']

if options.chroms:
    cmd += ['--chroms'] + options.chroms

cmd += ['--xg_paths', options.xg_basename + '.xg']
for name in options.names:
    cmd += [options.xg_basename + '_' + name + '.xg']

cmd += ['--gams']
for name in [options.xg_basename] + [options.xg_basename + '_' + n for n in options.names]:
    cmd += [os.path.join(options.out_store, 'aligned-{}-pe_default.gam'.format(os.path.basename(name)))]
cmd += ['--gam_names']
for name in [options.xg_basename] + [options.xg_basename + '_' + n for n in options.names]:
    cmd += [os.path.basename(name)]

cmd += ['--bams', os.path.join(options.outstore, 'bwa-mem-pe.bam')]
cmd += ['--bam_names', 'bwa-mem']

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
       '/venv/bin/aws', 's3', 'cp', log_name, 's3://{}'.format(os_log_name)]
print ' '.join(cmd)
subprocess.check_call(cmd)
