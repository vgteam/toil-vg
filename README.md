# TOIl-VG
## University of California, Santa Cruz Genomics Institute
### Please contact us on [github with any issues](https://github.com/BD2KGenomics/toil-vg/issues/new)

[vg](https://github.com/vgteam/vg) is a toolkit for DNA sequence analysis using varition graphs.  Toil-vg is a [toil](https://github.com/BD2KGenomics/toil)-based framework for running common vg pipelines at scale, either locally or on a distributed computing environment: 

`toil-vg run`: Given input graphs (one per chromosome) and reads (fastq file), produce a graph index (index can also be input), graph alignment (GAM), and VCF variant calls. 

`toil-vg index`: Produce an index from input graph(s).

`toil-vg map`: Produce graph alignment (gam) from input reads and index

`toil-vg call`: Produce VCF from input index and alignement (for single chromosome)

## Basic start-to-finish run on a local mac or linux machine

### Install toil-vg in a Python virtual environment

    virtualenv toilvenv
    source toilvenv/bin/activate
    git clone https://github.com/BD2KGenomics/toil-vg.git
    pip install --process-dependency-links ./toil-vg/
    pip install boto

### Install and setup docker on local machine and run the docker daemon
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
* It can be useful, especially for developers to run commands directly without docker.  This can be done by useing the --no_docker flag when running toil-vg, but all command-line tools (vg, bcftools, samtools, etc) must all be runnable from the command line

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
    toil-vg run ./local-toilvg-jobstore ./test_vg_input/NA12877.brca1.brca2.bam.fq NA12877 ./local-toilvg-output --graphs ./test_vg_input/BRCA1_BRCA2_may6.vg --chroms 13 --index_cores 2 --alignment_cores 2 --calling_cores 2 


## Basic start-to-finish run of the Toi-VG pipeline on Amazon EC2

### Install and setup local cgcloud

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
     

### Startup a Toil EC2 cluster with 4 m4.xlarge nodes

    cgcloud create-cluster toil -s 3 -t m4.xlarge --cluster-name toil-setup-test
    cgcloud ssh --admin -c toil-setup-test toil-leader 'sudo apt-get install -y git python-virtualenv && sudo pip install awscli'
    cgcloud ssh-cluster --admin --cluster-name toil-setup-test toil 'sudo pip install toil[aws,mesos]==3.5.0a1.dev251 boto3'

### Log onto leader node and set up dependencies

    cgcloud ssh -c toil-setup-test toil-leader
    virtualenv --system-site-packages toilvenv
    source toilvenv/bin/activate
    git clone --recursive https://github.com/BD2KGenomics/toil-vg.git
    pip install --process-dependency-links /home/mesosbox/toil-vg/

### Get test input files from s3
Any vg greaph (with primary reference path) and corresponding reads can be used here.   

    mkdir test_input/
    aws s3 cp s3://cgl-pipeline-inputs/vg_cgl/ci/BRCA1_BRCA2_may6.vg ./test_vg_input/BRCA1_BRCA2_may6.vg
    aws s3 cp s3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.brca1.brca2.bam.fq ./test_vg__input/NA12877.brca1.brca2.bam.fq

### Example run of VG toil-pipeline for variant calling on both chromosome 13 and 17 for sample NA12877

    toil-vg run aws:us-west-2:cmarkello-s3-jobstore ./test_vg_input/NA12877.brca1.brca2.bam.fq NA12877 ./local-toilvg-output --graphs ./test_vg_input/BRCA1_BRCA2_may6.vg --chroms 13 17 --index_cores 2 --alignment_cores 2 --calling_cores 2 --realTimeLogging --logInfo --batchSystem mesos --mesosMaster=mesos-master:5050
    toil clean aws:us-west-2:cmarkello-s3-jobstore 
    
### Run unit tests
(Requires credentials to download from S3 buckets)

From toil-vg directory, run `make test`
    
#### Delete cgcloud cluster

    cgcloud terminate-cluster --cluster-name toil-setup-test toil

## Running on Microsoft Azure

### Startup a Toil Azure cluster

- Go to the toil azure readme instructions [here](https://github.com/BD2KGenomics/toil/tree/master/contrib/azure#mesos-cluster-with-toil).
- Click on the `Deploy to Azure` button to use the toil Azure template for setting up a toil Azure cluster.
- Follow the instructions on setup. You will need to specify the following parameters as defined [here](https://github.com/BD2KGenomics/toil/tree/master/contrib/azure#template-parameters).
- Login to the master node of the Azure cluster by running `ssh <your_user_name>@<your_cluster_name>.<your_zone>.cloudapp.azure.com -p 2211`.

Remaing steps are identical to AWS, beginning with **Log onto leader node and set up dependencies** section above
