# TOIl-VG
## University of California, Santa Cruz Genomics Institute
### Please contact us on [github with any issues](https://github.com/BD2KGenomics/toil-vg/issues/new)

[vg](https://github.com/vgteam/vg) is a toolkit for DNA sequence analysis using varition graphs.  Toil-vg is a [toil](https://github.com/BD2KGenomics/toil)-based framework for running common vg pipelines at scale, either locally or on a distributed computing environment: 

`toil-vg run`: Given input graphs (one per chromosome) and reads (fastq file), produce a graph index (index can also be input), graph alignment (GAM), VCF variant calls, and (optionally) VCF comparison results. 

`toil-vg index`: Produce an index from input graph(s).

`toil-vg map`: Produce graph alignment (gam) from input reads and index

`toil-vg call`: Produce VCF from input index and alignement (for single chromosome)

## Installation

### Pip Installation

Installation requires Python.  We recommend installing within vituralenv as follows

    virtualenv toilvenv
    source toilvenv/bin/activate
    pip install toil-vg

### Docker
#### On linux
* Go to the main docker site and follow the instructions for the relevant linux distribution [here](https://docs.docker.com/engine/installation/linux/)
* Test to see if the docker daemon is running by running `docker version`
* If running `docker version` doesn't work, try adding `user` to docker group.
    * `sudo usermod -aG docker $USER`
    * log out and log back in

#### On Mac
* Install docker via instrictions found [here](https://docs.docker.com/docker-for-mac/)
* Test to see if the docker daemon is running by running `docker version`
* If running `docker version` doesn't work, try adding docker environment variables
    * `docker-machine start`
    * `docker-machine env`
    * `eval "$(docker-machine env default)"`
    
#### Running without Docker
* It can be useful, especially for developers to run commands directly without docker.  This can be done by useing the `--no_docker` flag when running toil-vg, but all command-line tools (vg, bcftools, samtools, etc) must all be runnable from the command line.  


## Configuration

A configuration file can be used as an alternative to most command line options.  A default configuration file can be generated using

    toil-vg generate-config > config.yaml

Pass this file to `toil-vg` commands using the `--config` option.

For non-trivial inputs, care must be taken to specify the resource requirements for the different pipeline phases (via the command line or by editing the config file), as they all default to single-core and 4G of ram.

    toil-vg generate-config --whole_genome > config_wg.yaml

## Testing

    make test

A faster test (but still several minutes) to see if toil-vg runs on the current machine.  Replace myname with a unique prefix: 

    ./bakeoff.sh -f myname f1.tsv

Or on a Mesos cluster

    ./bakeoff.sh -fm myname f1.tsv

In both cases, verify that f1.tsv contains a number (should be approx. 0.9).  Note that this script will create some directories (or S3 buckets) of the form `myname-bakeoff-out-store-brca1` and `myname-bakeoff-job-store-brca1`.  These will have to be manually removed. 

## A Note on IO conventions

The jobStore and outStore arguments to toil-vg are directories that will be created if they do not already exist.  When starting a new job, toil will complain if the jobStore exists, so use `toil clean <jobStore>` first.  When running on Mesos, these stores should be S3 buckets.  They are specified using the following format aws:region:bucket (see examples below).

All other input files can either either be local (best to specifiy absolute path) or URLs specified in the normal manner, ex : htpp://address/input_file or s3://bucket/input_file.  The config file must always be local.


## Running on Amazon EC2 with cgcloud

### Install and setup cgcloud

For more information on the cgcloud core tools, you can find them [here](https://github.com/BD2KGenomics/cgcloud/blob/master/README.md).
For more information on the cgcloud plugin for Toil and setting up AWS credentials, you can read about it [here](https://github.com/BD2KGenomics/cgcloud/blob/master/toil/README.rst).
For more information on the latest release of Toil, you can find the documentation [here](http://toil.readthedocs.io/en/latest/).

    sudo apt-get install python-virtualenv
    virtualenv ~/cgcloud

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

Create a test cluster with 2 m4.2xlarge nodes, and install toil-vg. You can modify the 2nd and 3rd lines below to install additional software on the leader node and whole cluster, respectively.  You will have to type `yes` when prompted for connecting to each node (the first time) to add it to known_hosts. 

    cgcloud create-cluster toil -s 1 -t m4.2xlarge --cluster-name toil-setup-test
    cgcloud ssh --admin -c toil-setup-test toil-leader 'sudo apt-get install -y git'
    cgcloud ssh-cluster --admin --cluster-name toil-setup-test toil 'sudo pip install toil-vg'

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
    cgcloud ssh-cluster --admin --cluster-name toil-index-cluster toil 'sudo pip install toil-vg'

Log on and switch to large disk volume.  It is best to run jobs within screen. 

    cgcloud ssh -c toil-index-cluster toil-leader
    screen
    cd /mnt/ephemeral/var/lib/mesos/
    toil-vg generate-config --whole_genome > wg.yaml

If they aren't available via a URL, download the input (chopped, common id space) .vg graphs, as [created here for example](https://github.com/vgteam/vg/wiki/working-with-a-whole-genome-variation-graph): 

Run the indexing (will take a couple days).  Make sure to edit the jobstore and output store agurments. Note that this invocation assumes 23 chromosome vg graphs are present in the current directory with names 1.vg, 2.vg ... X.vg, Y.vg.  Edit the `--graphs` and `--chroms` arguments to change. 

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

    cglcoud create-cluster toil -s 8 --instance-type r3.8xlarge  --leader-instance-type r3.2xlarge --cluster-name toil-map-cluster
    cgcloud ssh --admin -c toil-mac-cluster toil-leader 'sudo apt-get install -y aria2'    
    cgcloud ssh-cluster --admin --cluster-name toil-map-cluster toil 'sudo pip install toil-vg'

Log on and switch to large disk volume

    cgcloud ssh -c toil-map-cluster toil-leader
    screen
    cd /mnt/ephemeral/var/lib/mesos/
    toil-vg generate-config --whole_genome > wg.yaml

If they aren't available via a URL, download the input reads fastq (or fastq.gz) file using `aria2c -s 10 -x 10`.  If it is paired end, add `--interleaved` to the command below.  

Run the mapping.  Make sure to edit the jobstore and output store agurments, as well as the input index arguments to reflect the correct locations

    toil-vg run aws:us-west-2:myname-s3-jobstore ./reads.fastq.gz SAMPLE_NAME aws:us-west-2:myname-s3-outstore --workDir /mnt/ephemeral/var/lib/mesos/  --batchSystem=mesos --mesosMaster=mesos-master:5050 --gcsa_index s3://my-s3-outstore/my_index.gcsa --xg_index s3://my-s3-outstore/my_index.xg --id_ranges s3://my-s3-outstore/my_index_id_ranges.tsv  --realTimeLogging --logInfo --config wg.yaml --index_name my_index 2> index.log

If successful, this will produce gam files for each chromsome, as well as one VCF in the S3 output store

Terminate the cluster

    cgcloud terminate-cluster --cluster-name toil-map-cluster toil

## Running on Microsoft Azure

### Startup a Toil Azure cluster

- Go to the toil azure readme instructions [here](https://github.com/BD2KGenomics/toil/tree/master/contrib/azure#mesos-cluster-with-toil).
- Click on the `Deploy to Azure` button to use the toil Azure template for setting up a toil Azure cluster.
- Follow the instructions on setup. You will need to specify the following parameters as defined [here](https://github.com/BD2KGenomics/toil/tree/master/contrib/azure#template-parameters).
- Login to the master node of the Azure cluster by running `ssh <your_user_name>@<your_cluster_name>.<your_zone>.cloudapp.azure.com -p 2211`.

Remaing steps are identical to AWS, beginning with **Log onto leader node and set up dependencies** section above

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
