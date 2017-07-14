# TOIl-VG
## University of California, Santa Cruz Genomics Institute
### Please contact us on [github with any issues](https://github.com/BD2KGenomics/toil-vg/issues/new)

[vg](https://github.com/vgteam/vg) is a toolkit for DNA sequence analysis using variation graphs.  Toil-vg is a [toil](https://github.com/BD2KGenomics/toil)-based framework for running common vg pipelines at scale, either locally or on a distributed computing environment: 

`toil-vg run`: Given input vg graph(s), create indexes, map reads, then produce VCF variant calls.

`toil-vg index`: Produce a GCSA and/or XG index from input graph(s).

`toil-vg map`: Produce a graph alignment (GAM) for each chromosome from input reads and index

`toil-vg call`: Produce VCF from input XG index and GAM(s).

## Installation

### TOIL-VG Pip Installation

Installation requires Python and Toil.  We recommend installing within virtualenv as follows

    virtualenv --system-site-packages toilvenv
    source toilvenv/bin/activate
    pip install -I -U 'toil[aws,mesos]' toil-vg

### Docker

toil-vg can run vg, along with some other tools, via [Docker](http://www.docker.com).  Docker can be installed locally (not required when running via cgcloud), as follows. 
* [**Linux Docker Installation**](https://docs.docker.com/engine/installation/linux/): If running `docker version` doesn't work, try adding user to docker group with `sudo usermod -aG docker $USER`, then log out and back in.
* [**Mac Docker Installation**](https://docs.docker.com/docker-for-mac/): If running `docker version` doesn't work, try adding docker environment variables: `docker-machine start; docker-machine env; eval "$(docker-machine env default)"`
* **Running Without Docker**: If Docker is not installed or is disabled with `--container None`, toil-vg requires the following command line tools to be installed on the system: `vg, pigz, bcftools, tabix`.  `jq, samtools and rtg vcfeval` are also necessary for certain tests. 
    

## Configuration

A configuration file can be used as an alternative to most command line options.  A default configuration file can be generated using

    toil-vg generate-config > config.yaml

Pass this file to `toil-vg` commands using the `--config` option.

For non-trivial inputs, care must be taken to specify the resource requirements for the different pipeline phases (via the command line or by editing the config file), as they all default to single-core and 4G of ram.

To generate a default configuration for running at genome scale on a cluster with 32-core worker nodes, use

    toil-vg generate-config --whole_genome > config_wg.yaml

## Testing

    make test

A faster test to see if toil-vg runs on the current machine (Replace myname with a unique prefix): 

    ./bakeoff.sh -f myname f1.tsv

Or on a Toil cluster

    ./bakeoff.sh -fm myname f1.tsv

In both cases, verify that f1.tsv contains a number (should be approx. 0.9).  Note that this script will create some directories (or S3 buckets) of the form `myname-bakeoff-out-store-brca1` and `myname-bakeoff-job-store-brca1`.  These will have to be manually removed. 

## A Note on IO conventions

The jobStore and outStore arguments to toil-vg are directories that will be created if they do not already exist.  When starting a new job, toil will complain if the jobStore exists, so use `toil clean <jobStore>` first.  When running on Mesos, these stores should be S3 buckets.  They are specified using the following format aws:region:bucket (see examples below).

All other input files can either either be local (best to specify absolute path) or URLs specified in the normal manner, ex : http://address/input_file or s3://bucket/input_file.  The config file must always be local.  When using an S3 jobstore, it is preferable to pass input files from S3 as well, as they load much faster and less cluster time will be wasted importing data. 

## Running on Amazon EC2 with Toil

### Install Toil

Please read Toil's [installation documentation](http://toil.readthedocs.io/en/latest/install/basic.html)

Install the latest Toil prerelease locally.  This can be done with virtualenv as follows:

    virtualenv ~/toilvenv
    . ~/toilvenv/bin/activate
    pip install --pre toil[aws,mesos]

### Create a leader node

Please read Toil's cloud documentation [here](http://toil.readthedocs.io/en/latest/install/cloud.html) and [here](http://toil.readthedocs.io/en/latest/running/cloud.html)

Set your region:

	 export TOIL_AWS_ZONE="us-west-2a"	

Create an instance for the toil leader node.  For UCSC CGL members, the keypair name is typically the email address you use to log into AWS:

    toil launch-cluster myleader --nodeType=t2.micro --keyPairName <your AWS keypair name>

Log into the leader:

    toil ssh-cluster myleader

In order to use `screen` (recommended for long jobs), you need to run `script` first.

In order to log onto a worker node instead of the leader, find its public IP from the EC2 Management Console or command line, and log in using the core username: `ssh core@public-ip`

Install toil-vg on the leader following [Toil's Hot-deployment instructions](http://toil.readthedocs.io/en/latest/deploying.html#hot-deploying-toil)

    mkdir work
    cd work
    virtualenv --system-site-packages venv
    . venv/bin/activate
    pip install toil-vg

Destroy the leader when finished with it.  After logging out with `exit`:

    toil destroy-cluster myleader

### Small AWS Test

Run a small test from the leader node as follows.  

    wget https://raw.githubusercontent.com/BD2KGenomics/toil-vg/master/bakeoff.sh
    chmod u+x ./bakeoff.sh
    ./bakeoff.sh -fm <NAME>

### Processing a Whole Genome 

From the leader node, begin by making a toil-vg configuration file suitable for processing whole-genomes, then customizing it as necessary.

    toil-vg generate-config --whole_genome > wg.yaml

Assuming a series of input vg graphs, one per chromosome, as [created here for example](https://github.com/vgteam/vg/wiki/working-with-a-whole-genome-variation-graph) exists on GRAPH_LOCATION (can be any valid Toil path), files will be written to the S3 bucket, OUT_STORE and the S3 bucket, JOB_STORE, will be used by Toil (both buckets created automatically if necessary; do not prefix OUT_STORE or JOB_STORE with s3://):

    MASTER_IP=`ifconfig eth0 |grep "inet addr" |awk '{print $2}' |awk -F: '{print $2}'`
    toil-vg index aws:us-west-2:JOB_STORE aws:us-west-2:OUT_STORE --workDir /var/lib/toil/  --batchSystem=mesos --mesosMaster=${MASTER_IP}:5050  --graphs $(for i in $(seq 22; echo X; echo Y); do echo GRAPH_LOCATION/${i}; done) --chroms $(for i in $(seq 22; echo X; echo Y); do echo $i; done) --realTimeLogging --logInfo --config wg.yaml --index_name my_index --defaultPreemptable --preemptableNodeType i3.8xlarge:1.00 --maxPreemptableNodes 5 --nodeType i3.8xlarge --provisioner aws 2> index.log

Note that the spot request node type (i3.8xlarge) and amount ($1.00) can be adjusted in the above command.  Keep in mind that indexing is very memory and disk intensive.

If successful, this will produce for files in s3://OUT_STORE/

    my_index.xg
    my_index.gcsa
    my_index.gcsa.lcp
    my_index_id_ranges.tsv

We can now align reads and produce a VCF in a single call to `toil-vg run`. (see `toil-vg map` and `toil-vg call` to do separately).  The invocation is similar to the above, except we use r3.8xlarge instances as we do not need as much disk and memory.

    toil-vg run aws:us-west-2:JOB_STORE READ_LOCATION/reads.fastq.gz SAMPLE_NAME aws:us-west-2:OUT_STORE --workDir /var/lib/toil/  --batchSystem=mesos --mesosMaster=${MASTER_IP}:5050 --gcsa_index s3://OUT_STORE/my_index.gcsa --xg_index s3://OUT_STORE/my_index.xg --id_ranges s3://${OUT_STORE}/my_index_id_ranges.tsv  --realTimeLogging --logInfo --config wg.yaml --index_name my_index --interleaved --defaultPreemptable --preemptableNodeType r3.8xlarge:0.85 --maxPreemptableNodes 25 --nodeType r3.8xlarge --provisioner aws 2> map_call.log

If successful, this command will create a VCF file as well as a GAM for each input chromosome in s3://OUT_STORE/


## Running on Amazon EC2 with cgcloud

**Note: Unfortunately, cgcloud is no longer being maintained.  Many of its images are becoming obselete or vanishing entirely.  This section is deprecated until further notice.  toil-vg can now be run on AWS directly from Toil using hte latter's built-in AWS provisioner (see above).**

### Install and setup cgcloud

Please see [here for more information on the cgcloud core tools](https://github.com/BD2KGenomics/cgcloud/blob/master/README.md), [here for more information on the cgcloud plugin for Toil and setting up AWS credentials](https://github.com/BD2KGenomics/cgcloud/blob/master/toil/README.rst), and [here for more information on the latest release of Toil](http://toil.readthedocs.io/en/latest/).

    sudo apt-get install python-virtualenv
    virtualenv ~/cgcloud
	 pip install cgcloud-toil

Edit your ~/.profile or ~/.bash_profile and add the following lines:

    export CGCLOUD_ZONE=us-west-2a`
    export CGCLOUD_PLUGINS="cgcloud.toil:$CGCLOUD_PLUGINS"

Setup credentials for your AWS account in `~/.aws/credentials`:
    
    [default]
    aws_access_key_id=PASTE_YOUR_FOO_ACCESS_KEY_ID_HERE
    aws_secret_access_key=PASTE_YOUR_FOO_SECRET_KEY_ID_HERE
    region=us-west-2
    source ~/cgcloud/bin/activate
    pip install cgcloud-toil
    cgcloud register-key ~/.ssh/id_rsa.pub
 
Create a cluster image:
 
    cgcloud create -IT toil-box

### Test cgcloud

Create a test cluster with 2 m4.2xlarge nodes, and install toil-vg. You can modify the 2nd and 3rd lines below to install additional software on the leader node and whole cluster, respectively.  You will have to type `yes` when prompted for connecting to each node (the first time) to add it to known_hosts. We reinstall toil on the cluster to make sure all versions are up to date and consistent (cgcloud images currently contain outdated version of Toil)

    cgcloud create-cluster toil -s 1 -t m4.2xlarge --cluster-name toil-setup-test
    cgcloud ssh --admin -c toil-setup-test toil-leader 'sudo apt-get install -y git'
    cgcloud ssh-cluster --admin --cluster-name toil-setup-test toil "sudo pip install -I 'toil[aws,mesos]>=3.6.0' toil-vg"

Login to the leader node and run a test:

    cgcloud ssh -c toil-setup-test toil-leader
    git clone --recursive https://github.com/BD2KGenomics/toil-vg.git
    toil-vg/bakeoff.sh -fm myname output.tsv
    cat output.tsv

Terminate the cluster:

    exit (if still logged into leader node)
    cgcloud terminate-cluster --cluster-name toil-setup-test toil

Note that bakeoff.sh script will create some S3 buckets of the form `myname-bakeoff-out-store-brca1` and `myname-bakeoff-job-store-brca1`.  These will have to be manually removed. 

### Run a whole genome on AWS

#### Indexing

Indexing requires lots of storage and RAM, and much of it cannot be distributed (single call to `vg index`).  We therefore use a single node.

    cgcloud create-cluster toil -s 1 --instance-type i2.8xlarge  --leader-instance-type r3.xlarge --cluster-name toil-index-cluster
    cgcloud ssh-cluster --admin --cluster-name toil-index-cluster toil "sudo pip install -I 'toil[aws,mesos]>=3.6.0' toil-vg"

Log on and switch to large disk volume.  It is best to run jobs within screen. 

    cgcloud ssh -c toil-index-cluster toil-leader
    screen
    cd /mnt/ephemeral/var/lib/mesos/
    toil-vg generate-config --whole_genome > wg.yaml

If they aren't available via a URL, download the input (chopped, common id space) .vg graphs, as [created here for example](https://github.com/vgteam/vg/wiki/working-with-a-whole-genome-variation-graph): 


Run the indexing (will take about 40 hours).  **Make sure to edit the jobstore and output store arguments to change "myname"**. Note that this invocation assumes 24 chromosome vg graphs are present in the current directory with names 1.vg, 2.vg ... X.vg, Y.vg.  Edit the `--graphs` and `--chroms` arguments to change. 

    toil-vg index aws:us-west-2:myname-s3-jobstore aws:us-west-2:myname-s3-outstore --workDir /mnt/ephemeral/var/lib/mesos/  --batchSystem=mesos --mesosMaster=mesos-master:5050  --graphs $(for i in $(seq 22; echo X; echo Y); do echo /mnt/ephemeral/var/lib/mesos/$i.vg; done) --chroms $(for i in $(seq 22; echo X; echo Y); do echo $i; done) --realTimeLogging --logInfo --config wg.yaml --index_name my_index 2> index.log

If successful, this will produce for files in the S3 output store

    my_index.xg
    my_index.gcsa
    my_index.gcsa.lcp
    my_index_id_ranges.tsv

Terminate the cluster

    cgcloud terminate-cluster --cluster-name toil-index-cluster toil
    
#### Mapping and variant calling

This part of the pipeline is more distributed, so we make a larger cluster with less storage.  Note: indexing, mapping, and calling can be done in a single invocation of `toil-vg run` (by leaving out `--gcsa_index`, `--xg_index`, and `--id_ranges`) if you feel your setup is up to it.   

    cgcloud create-cluster toil -s 8 --instance-type r3.8xlarge  --leader-instance-type r3.2xlarge --cluster-name toil-map-cluster
    cgcloud ssh --admin -c toil-map-cluster toil-leader 'sudo apt-get install -y aria2'    
    cgcloud ssh-cluster --admin --cluster-name toil-map-cluster toil "sudo pip install -I 'toil[aws,mesos]>=3.6.0' toil-vg"

Log on and switch to large disk volume

    cgcloud ssh -c toil-map-cluster toil-leader
    screen
    cd /mnt/ephemeral/var/lib/mesos/
    toil-vg generate-config --whole_genome > wg.yaml

It is best to pass in large input files via S3 when possible, but if they aren't available via a URL, download the input reads fastq (or fastq.gz) file using `aria2c -s 10 -x 10 url`.  If the reads are not paired-end, remove `--interleaved` from the command below.  

Run the mapping.  **Make sure to edit the jobstore and output store agurments, as well as the input index and reads arguments and to reflect the correct locations**

    toil-vg run aws:us-west-2:myname-s3-jobstore ./reads.fastq.gz SAMPLE_NAME aws:us-west-2:myname-s3-outstore --workDir /mnt/ephemeral/var/lib/mesos/  --batchSystem=mesos --mesosMaster=mesos-master:5050 --gcsa_index s3://myname-s3-outstore/my_index.gcsa --xg_index s3://myname-s3-outstore/my_index.xg --id_ranges s3://myname-s3-outstore/my_index_id_ranges.tsv  --realTimeLogging --logInfo --config wg.yaml --index_name my_index --interleaved 2> map.log

If successful, this will produce a gam file for each chromsome, as well as a whole-genome VCF in the S3 output store


Terminate the cluster

    cgcloud terminate-cluster --cluster-name toil-map-cluster toil

## Running on Microsoft Azure

### Startup a Toil Azure cluster

- Go to the toil azure readme instructions [here](https://github.com/BD2KGenomics/toil/tree/master/contrib/azure#mesos-cluster-with-toil).
- Click on the `Deploy to Azure` button to use the toil Azure template for setting up a toil Azure cluster.
- Follow the instructions on setup. You will need to specify the following parameters as defined [here](https://github.com/BD2KGenomics/toil/tree/master/contrib/azure#template-parameters).
- Login to the master node of the Azure cluster by running `ssh <your_user_name>@<your_cluster_name>.<your_zone>.cloudapp.azure.com -p 2211`.

The remaining steps are identical to running on AWS

## Local test without AWS data

### Get test input files from Synapse
* Register for a [Synapse account](https://www.synapse.org/#!RegisterAccount:0)
* Create directory to store input files: `mkdir ./test_vg_input/`
* Either download the samples from the [website GUI](https://www.synapse.org/#!Synapse:syn7562100) or use the Python API to download into the `./test_vg_input/` directory
* `pip install synapseclient`
* `python`
    * `import synapseclient`
    * `syn.login('foo@bar.com', 'password')`
    * Get the test vg graph
        * `syn.get('syn7562120', downloadLocation='./test_vg_input/')`
    * Get the test read set
        * `syn.get('syn7562136', downloadLocation='./test_vg_input/')`

### Example run for BRCA2 alignment and variant calling for sample NA12877
    rm -fr ./local-toilvg-output
    toil clean ./local-toilvg-jobstore
    toil-vg run ./local-toilvg-jobstore ./test_vg_input/NA12877.brca1.brca2.bam.fq NA12877 ./local-toilvg-output --graphs ./test_vg_input/BRCA1_BRCA2_may6.vg --chroms 13 
