# TOIl-VG
[UCSC Computational Genomics Lab](https://cgl.genomics.ucsc.edu/)

[vg](https://github.com/vgteam/vg) is a toolkit for DNA sequence analysis using variation graphs.  Toil-vg is a [Toil](http://toil.ucsc-cgl.org/)-based framework for running common vg pipelines at scale, either locally or on a distributed computing environment: 

* `toil-vg construct`: Create vg graphs and indexes from FASTA and VCF, constructing contigs in parallel.  
* `toil-vg index`: Produce a GCSA, GBWT and/or XG index from input vg graphs.
* `toil-vg map`: Produce a graph alignment (GAM) for each chromosome from input reads and index
* `toil-vg call`: Produce VCF from input XG index and GAM(s).

Why use toil-vg?

* Higher-level CLI simplifies tasks that require several complex bash and vg commands to do.  For example, constructing and indexing a graph from a 1000 Genomes VCF takes [dozens](https://github.com/vgteam/vg/wiki/Working-with-a-whole-genome-variation-graph) of [commands](https://github.com/jltsiren/gbwt/wiki/Construction-Benchmarks) using vg directly.  They can all (including those for downloading input data), and more, be replaced by a single call to `toil-vg construct`.
* In practice, it is necessary to distribute mapping and calling jobs across multiple compute nodes when working with human-sized genomes.  toil-vg takes care of this. 
* [Toil's](http://toil.ucsc-cgl.org/) support for running on the cloud, resuming lost or failed jobs and running vg and other dependencies via Docker are further reasons to run toil-vg on large datasets.
* toil-vg provides benchmarking scripts to run large vg experiments comparing different graphs, aligners, callers, mapping parameters, etc.

#### Please contact us here on [github with any issues](https://github.com/BD2KGenomics/toil-vg/issues/new)

#### See the [Wiki](https://github.com/vgteam/toil-vg/wiki) in addition to below for examples.

## Installation

Installation requires Python and [Toil](https://toil.readthedocs.io/en/latest/gettingStarted/install.html).  We recommend installing within a virtualenv as follows

    virtualenv toilvenv
    source toilvenv/bin/activate
    pip install toil[aws,mesos]==3.18.0
    pip install toil-vg
    
Developers may want to work with the latest master branch from github instead:

    git clone https://github.com/vgteam/toil-vg.git
    cd toil-vg
    virtualenv toilvenv
    source toilvenv/bin/activate
    make prepare
    make develop
    
    # After installing in this fashion, the tests can be run with
    make test
    
### Docker

toil-vg can run vg, along with some other tools, via [Docker](http://www.docker.com).  Docker can be installed locally as follows:
* [**Linux Docker Installation**](https://docs.docker.com/engine/installation/linux/): If running `docker version` doesn't work, try adding user to docker group with `sudo usermod -aG docker $USER`, then log out and back in.
* [**Mac Docker Installation**](https://docs.docker.com/docker-for-mac/): If running `docker version` doesn't work, try adding docker environment variables: `docker-machine start; docker-machine env; eval "$(docker-machine env default)"` **Note** If you get an error message of the form `APIError: 502 Server Error: Bad Gateway (Mounts denied:...`, you can fix this by either adding `/var/folders` to the directories list in the "File Sharing" section of the Docker preferences dialog, or by specifying a temporary working directory that already is shared using the `--tempDir` toil-vg option or the `TMPDIR` environment variable. 
* **Running Without Docker**: If Docker is not installed or is disabled with `--container None`, toil-vg requires the following command line tools to be installed on the system: `vg, pigz, bcftools, tabix`.  `jq`, `samtools`, `rtg vcfeval`, 'hap.py', and 'Platypus.py' are also necessary for certain tests. 
    

## Configuration

A configuration file can be used as an alternative to most command line options.  A default configuration file can be generated using

    toil-vg generate-config > config.yaml

Pass this file to `toil-vg` commands using the `--config` option.

For non-trivial inputs, care must be taken to specify the resource requirements for the different pipeline phases (via the command line or by editing the config file), as they all default to single-core and 4G of ram.

To generate a default configuration for running at genome scale on a cluster with 32-core worker nodes, use

    toil-vg generate-config --whole_genome > config_wg.yaml

## A Note on IO conventions

The jobStore and outStore arguments to toil-vg are directories that will be created if they do not already exist.  When starting a new job, toil will complain if the jobStore exists, so use `toil clean <jobStore>` first.  When running on Mesos, these stores should be S3 buckets.  They are specified using the following format aws:region:bucket (see examples below).

All other input files can either either be local (best to specify absolute path) or URLs specified in the normal manner, such as: `http://address/input_file` or `s3://bucket/input_file`.  The config file must always be local.  When using an S3 jobstore, it is preferable to pass input files from S3 as well, as they load much faster and less cluster time will be wasted importing data. 

## Example: Construct and index 1000 Genomes HS38D1 graph

The following command will construct and index a graph from the 1000 Genomes calls for the HS38D1 (GRCH38 + decoys) reference.  

```
# Toil boilerplate
export TOIL_OPTS="--realTimeLogging --logInfo --realTimeStderr"
# Directory for Toil's temporary files
export TOIL_JS="./my-jobstore"
# All output will be written here
export TOIL_OS="./my-output"

# Construct graph and all (XG, GBWT, GCSA, Snarls, id-ranges) indexes

toil-vg construct $TOIL_JS $TOIL_OS --pangenome --out_name 1kg_hs38d1 --all_index --merge_graphs --fasta_regions --add_chr_prefix --whole_genome_config $TOIL_OPTS --fasta ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa --vcf $(for i in $(seq 1 22; echo X; echo Y); do echo ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${i}_GRCh38.genotypes.20170504.vcf.gz; done) --regions $(for i in $(seq 1 22; echo X; echo Y); do echo chr${i}; done)
```

Filtering out low-frequency variants (often recommended) can be done by replacing `--pangenome` by `--min_af 0.01` (both options can be used to simultaneously create filtered and unfiltered graphs).  This command uses `--fasta_regions` to add every sequence from the input fasta file to the graph.  The `--regions_regex` option can be used to further fine-tune this.  For example, `--regions_regex 'chr[M,EBV]' 'chr.*decoy'` would only add chrM, chrEBV and decoys  (and not unplaced scaffolds).  It is good practice to use `--regions` to explicitly define the regions corresponding to each VCF in order, whether or not `--fasta_regions` is used.  

By default, this command will use every core on the system.  Use `--maxCores` to limit to fewer processes.  It is recommended to have 3TB disk, 256GB RAM and 32 cores to run this, and even then it will take several days to complete. See below for how the above command can be adapted to run on Amazon EC2.

## Running on Amazon EC2 with Toil

### Create a leader node

    wget https://raw.githubusercontent.com/BD2KGenomics/toil-vg/master/scripts/create-ec2-leader.sh
    ./create-ec2-leader.sh <leader-name> <keypair-name>

Log into the leader with

    toil ssh-cluster <leader-name> --zone us-west-2a

In order to log onto a worker node instead of the leader, find its public IP from the EC2 Management Console or command line, and log in using the core username: `ssh core@public-ip`

Destroy the leader when finished with it.  After logging out with `exit`:

    toil destroy-cluster -z us-west-2a myleader

### Run a job

Log onto the leader node with `toil ssh-cluster` as described above and open a `screen` session.  The same construction example as above can now be run after adapting the input environment variables as follows:

```
# Use maximum of 8 i3.8xlarge nodes with spot bit of $0.80
export NODE_OPTS="--nodeTypes i3.8xlarge:0.80,i3.8xlarge --maxNodes 8"
# Set Toil to use AWS autoscaling
export AWS_OPTS="--defaultPreemptable --batchSystem mesos --provisioner aws --retryCount 3 --metrics"
# Toil boilerplate
export TOIL_OPTS="$NODE_OPTS $AWS_OPTS --realTimeLogging --logInfo --realTimeStderr"
# S3 Bucket for Toil's temporary files
export TOIL_JS="aws:us-west-2:my-jobstore"
# All output will be written here
export TOIL_OS="aws:us-west-2/my-bucket/hs38d1-output"
```

For a performance dashboard, browse to `localhost:3000` on the computer from which you ran `toil ssh-cluster` (this is enabled by the `--metrics` option above)
