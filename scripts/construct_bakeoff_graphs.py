#!/usr/bin/python
"""
Locally regenerate all the bakeoff regions graphs and indexes that are
found here s3://vg-data/bakeoff/
The input fasta's and vcf's are expected to be there already.
Assumes you have authenticated S3 access configured. If not, the files are
mirrored to https://courtyard.gi.ucsc.edu/~anovak/vg-data/bakeoff/ 
"""

import os, sys, subprocess


region_to_bed_hg38 = {
    'BRCA1':('17', 43044293, 43125482),
    'BRCA2':('13', 32314860, 32399849),
    'SMA':('5', 69216818, 71614443),
    'MHC':('6', 28510119, 33480577)
}

def get_vcf_coords_hg38(region):
    r = region_to_bed_hg38[region]
    # change bed to 1-based inclusive
    return '{}:{}-{}'.format(r[0], r[1] + 1, r[2])

def get_vcf_path_hg38(region):
    return 's3://vg-data/bakeoff/1kg_hg38-{}.vcf.gz'.format(region)

def get_fasta_path_hg38(region):
    chrom = region_to_bed_hg38[region][0]
    return 's3://vg-data/bakeoff/chr{}.fa.gz'.format(chrom)

if len(sys.argv) not in [3,4]:
    print "Usage: {} jobstore outstore <config>".format(sys.argv[0])
    sys.exit(1)

job_store = sys.argv[1]
out_store = sys.argv[2]
config = sys.argv[3] if len(sys.argv) == 4 else None
config_opts = [] if not config else ['--config', config]

for region in ['BRCA1', 'BRCA2', 'SMA', 'MHC']:
    # make the graphs/indexes and a bunch of controls
    cmd = ['toil-vg', 'construct', job_store, out_store,
           '--vcf', get_vcf_path_hg38(region),
           '--fasta', get_fasta_path_hg38(region),
           '--regions', get_vcf_coords_hg38(region),
           '--out_name', 'snp1kg-{}'.format(region),
           '--alt_paths', '--realTimeLogging',
           '--control_sample', 'HG00096',
           '--min_af', '0.0335570469',
           '--primary', '--filter_samples', 'HG00096',
           '--xg_index', '--gcsa_index', '--gbwt_index'] + config_opts
    
    subprocess.check_call(cmd)

    # make the names consistent to what we've been using
    for os_file in os.listdir(out_store):
        prefix = 'snp1kg-{}'.format(region)
        if os_file.startswith(prefix):
            if os_file.endswith('.gcsa.lcp'):
                ext = '.gcsa.lcp'
                name = os_file[:-len(ext)]
            else:
                name, ext = os.path.splitext(os_file)
            new_name = 'snp1kg' + name[len(prefix):] + '-{}'.format(region) + ext
            if new_name.startswith('snp1kg_primary'):
                new_name = new_name[len('snp1kg_'):]
            elif new_name.startswith('snp1kg_minaf_0.0335570469'):
                new_name = 'snp1kg_threshold10' + new_name[len('snp1kg_minaf_0.0335570469'):]
            elif new_name.startswith('snp1kg_filter'):
                new_name = 'snp1kg_filter_HG00096' + new_name[len('snp1kg_filter'):]
            if os_file != new_name:
                cmd = ['mv', os.path.join(out_store, os_file), os.path.join(out_store, new_name)]
                subprocess.check_call(cmd)



