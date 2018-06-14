#!/usr/bin/python2.7
"""
Simulate from a baseline graph made with construct-baseline-hs37d5-ec2.py
Currently needs to be run as follows: scripts/sim-ec2.py
"""
# Doesn't really simply toil-vg sim interface much, but does take care of naming
import os, sys, subprocess, argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, 
        formatter_class=argparse.RawDescriptionHelpFormatter)    
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("basename", help="input file prefix.  will look for [basename]_thread_0.xg,"
                        "[basename]_thread_1.xg, [basename].xg")
    parser.add_argument("num_reads", type=int, help="number of read pairs")
    parser.add_argument("--leader", help="name of leader created with create-ec2-leader.sh")    
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--fastq", help="template fastq (default giab nist sample).  use \"None\" to turn off",
                        default="ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/"
                        "NIST_NA12878_HG001_HiSeq_300x/131219_D00360_005_BH814YADXX/"
                        "Project_RM8398/Sample_U5a/U5a_AGTCAA_L002_R1_007.fastq.gz")    
    parser.add_argument("--sim_opts", help="override simulation options")
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")
    args = args[1:]        
    return parser.parse_args(args)
options = parse_args(sys.argv)

if options.fastq.lower() == "none":
    options.fastq = None
    default_sim_opts = None
else:
    default_sim_opts = '-p 570 -v 65 -i 0.002 -I'
sim_opts = options.sim_opts if options.sim_opts else default_sim_opts
    
def print_num(n):
    """ there must be a package that does this better. oh well"""
    if n >= 1000000000:
        return "{}{}".format(round(float(n) / 1000000000, 2), "B")    
    elif n >= 1000000:
        return "{}{}".format(round(float(n) / 1000000, 2), "M")
    elif n >= 1000:
        return "{}{}".format(round(float(n) / 1000, 2), "K")
    else:
        return str(n)

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)

out_name = os.path.basename(options.basename) + '_sim_{}_'.format(print_num(options.num_reads))
if options.fastq:
    out_name += 'trained'
else:
    out_name += 'ekg'
if sim_opts:
    out_name += sim_opts.replace(' ', '')
out_name += '-s{}'.format(options.seed)
    
log_name = 'sim_{}.log'.format(out_name)
os_log_name = os.path.join(options.out_store[options.out_store.rfind(':')+1:], os.path.basename(log_name))

cmd = ['sim', options.job_store, 
       options.basename + '_thread_0.xg',
       options.basename + '_thread_1.xg',
       str(options.num_reads),       
       options.out_store,
       '--logFile', log_name,
       '--out_name', out_name,
       '--annotate_xg', options.basename + '.xg',       
       '--seed', str(options.seed),
       '--sim_chunks', str(int(max(1, options.num_reads / 500000))),
       '--gam', '--fastq_out']
if options.fastq:
    cmd += ['--fastq', options.fastq]
if sim_opts:
    cmd += ['--sim_opts', '\"{}\"'.format(sim_opts)]

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
           '/venv/bin/aws', 's3', 'cp', log_name, 's3://{}'.format(os_log_name)]
    print ' '.join(cmd)
    subprocess.check_call(cmd)

