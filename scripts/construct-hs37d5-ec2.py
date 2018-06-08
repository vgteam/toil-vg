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
    parser.add_argument("job_store")
    parser.add_argument("out_store")
    parser.add_argument("--leader", help="name of leader created with create-ec2-leader.sh."
                        "run on current machine if not specified")    
    parser.add_argument("--chroms", nargs='+', default=None,
                        help="chromosome(s) (default: everything in fasta)")
    parser.add_argument("--config", help="path of config on leader")
    parser.add_argument("--restart", action="store_true", help="resume toil workflow")
    parser.add_argument("--gbwt", action="store_true", help="make gbwt")
    parser.add_argument("--xg", action="store_true", help="make xg")
    parser.add_argument("--gcsa", action="store_true", help="make gcsa")
    parser.add_argument("--snarls", action="store_true", help="make snarls")
    parser.add_argument("--sample", help="sample name to use for control/haplo/sample graphs")
    parser.add_argument("--pos_control", action="store_true", help="make positive control graph")
    parser.add_argument("--neg_control", action="store_true", help="make negative control graph")
    parser.add_argument("--sample_graph", action="store_true", help="make sample graph")
    parser.add_argument("--haplo_graph", action="store_true", help="make haplo thread graphs")
    parser.add_argument("--primary", action="store_true", help="make primary graph")
    parser.add_argument("--minaf", type=float, nargs='+', help="make min allele filter graph")
    parser.add_argument("--alt_paths", action="store_true", help="force alt paths")
    parser.add_argument("--filter_ceph", action="store_true", help="filter private CEPH variants")    
    parser.add_argument("--node", help="toil node type(s) (i3.8xlarge:0.90,i3.8xlarge). can be comma-separated list", default="i3.8xlarge:0.90,i3.8xlarge")
    parser.add_argument("--max_node", help="Max nodes for each type (can be comman-separated list)",
                        default="8,4")
    parser.add_argument("--out_name", help="Prefix output files with this name")
    args = args[1:]        
    return parser.parse_args(args)
options = parse_args(sys.argv)

if options.pos_control or options.neg_control or options.sample_graph or options.haplo_graph:
    assert options.sample

def get_vcf_paths_hs37d5(phased, chroms = None, sample = None):
    """ Get list of paths for HS37D5 vcfs. If not phased, then a whole-genome vcf will be returned.  If phased
    and chromosomes provided, a list of phased vcfs are returned.  If sample is HG002, a single GIAB whole-genome phased
    vcf is returned no matter what.
    """
    if sample == 'HG002':
        out_vcfs = ['ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/AshkenazimTrio/HG002_NA24385_son/NISTv3.3.2/GRCh37/HG002_GRCh37_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-22_v.3.3.2_highconf_triophased.vcf.gz']
    else:
        if not phased:
            out_vcfs = ['ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz']
        else:
            assert chroms
            base = 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
            out_vcfs = []
            for chrom in chroms:
                try:
                    if int(chrom) in range(1, 23):
                        out_vcfs.append(os.path.join(base, 'ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'.format(chrom)))
                except:
                    pass
                if chrom == 'X':
                    out_vcfs.append(os.path.join(base, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz'))
                elif chrom == 'Y':
                    out_vcfs.append(os.path.join(base, 'ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz'))
            assert len(out_vcfs) == len(chroms)
    return out_vcfs

def get_fasta_path_hs37d5():
    return 'ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz'

if not options.job_store.startswith('aws:'):
    options.job_store = 'aws:us-west-2:{}'.format(options.job_store)
if not options.out_store.startswith('aws:'):
    options.out_store = 'aws:us-west-2:{}'.format(options.out_store)    

if options.out_name:
    out_name = options.out_name
else:
    out_name = 'snp1kg' if not options.chroms or len(options.chroms) > 3 else 'snp1kg_{}'.format('_'.join(options.chroms))
log_name = 'construct_{}.log'.format(out_name)
os_log_name = os.path.join(options.out_store[options.out_store.rfind(':')+1:], os.path.basename(log_name))

cmd = ['construct', options.job_store, options.out_store,
       '--fasta', get_fasta_path_hs37d5(),
       '--out_name', out_name,
       '--logFile', log_name]

# Note config file path is on the leader!!!!  Should fix to copy it over, but not sure how.
cmd += ['--config', options.config] if options.config else ['--whole_genome_config']

if options.sample != 'HG002':
    cmd += ['--pangenome']

# do we want a vcf with sample phasing in it?
is_phased = options.gbwt or options.alt_paths or options.filter_ceph or options.sample
if options.minaf:
    cmd += ['--min_af', ' '.join([str(minaf) for minaf in options.minaf])]

if options.primary:
    cmd += ['--primary']

if options.xg:
    cmd += ['--xg_index']

if options.gcsa:
    cmd += ['--gcsa_index']

if options.gbwt:
    cmd += ['--gbwt_index', '--gbwt_prune']
    
if options.gbwt or options.alt_paths or options.sample_graph or options.haplo_graph:
    cmd += ['--alt_paths']

if options.snarls:
    cmd += ['--snarls_index']

if options.pos_control:
    cmd += ['--pos_control', options.sample]

if options.neg_control:
    cmd += ['--neg_control', options.sample]

if options.sample_graph:
    cmd += ['--sample_graph', options.sample]

if options.haplo_graph:
    cmd += ['--haplo_sample', options.sample]

if options.filter_ceph:
    cmd += ['--filter_ceph']

if options.chroms:
    # restrict to specified chromosome(s)
    cmd += ['--regions'] + options.chroms
    vcfs = get_vcf_paths_hs37d5(is_phased, options.chroms, options.sample)
else:
    # do all chromsomes as well as decoys
    cmd += ['--fasta_regions']
    vcfs = get_vcf_paths_hs37d5(is_phased, [chrom for chrom in range(1, 23) + ['X', 'Y']], options.sample)
for vcf in vcfs:
    cmd += ['--vcf', vcf]

if options.restart:
    cmd += ['--restart']
else:
    subprocess.check_call(['toil', 'clean', options.job_store])


ec2_cmd = ['scripts/ec2-run.sh', '-n', options.node, '-m', options.max_node]
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

