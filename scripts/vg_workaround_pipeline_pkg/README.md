# NIH Biowulf VG Pipeline
## University of California, Santa Cruz Genomics Institute
### Please contact us on [github with any issues](https://github.com/BD2KGenomics/toil-vg/issues/new)
[vg](https://github.com/vgteam/vg) is a toolkit for DNA sequence analysis using variation graphs.  Toil-vg is a [toil](https://github.com/BD2KGenomics/toil)-based framework for running common vg pipelines at scale, either locally or on a distributed computing environment.

The pipeline runs the mapping and surjection of the VG graph alignment process intended for use on the NIH Biowulf cluster. It is mainly written to workaround the temporary incompatibility of the toil-based Toil-vg pipeline with the Slurm batch system that the NIH HPC systems use.
This package is only for use on NIH Biowulf systems. Further information on how to run various jobs on Biowulf instances can be found [here](https://hpc.nih.gov/).

## Basic Installation

First download the github repo
    
    git clone https://github.com/vgteam/toil-vg.git
    
Next add the bash scripts to your `PATH` environment variable

    export PATH=$PATH:toil-vg/scripts/vg_workaround_pipeline_pkg
    
## Basic WGS runs on an individual sample

To run an individual sample in this pipeline, you will first need to collect:
- Read files in fasta format (paired-end is only supported for now)
- A directory that contains VG graph reference [XG](https://github.com/vgteam/vg/wiki/File-Formats#xg-xg-lightweight-graph--path-index) and [GCSA](https://github.com/vgteam/vg/wiki/File-Formats#gcsa-gcsa-generalized-compressed-suffix-array-index) indexes.
- A work directory with plenty of disk space (~1TB for a single WGS sample).
- A tag for the version of VG which is used for importing the VG container that runs the VG program through this pipeline. Tag names can be aquired through the following [quay.io repo](https://quay.io/repository/vgteam/vg?tag=latest&tab=tags). An example tag from this repo would be "v1.6.0-939-g9fa94410-t149-run". It's preferable that you use a contain tag that ends in "run" instead of "build" for more efficient runs. 

To execute the pipeline, run the following command:

    sample_master_script.sh -i ${SAMPLE_NAME} -r ${INPUT_READ_FILE_1} -r ${INPUT_READ_FILE_2} -g ${GRAPH_FILES_DIR_PATH} -w ${WORK_DIR} -c ${VG_CONTAINER} > ${WORK_DIR}/master_script.stdout
    
### An example run of a wgs sample
    
    SAMPLE_NAME="NA12878
    READ_DIR="/data/markellocj/CEPH_TRIO_run/test_reads_${COHORT_NAME}"
    INPUT_READ_FILE_1="${READ_DIR}/${SAMPLE_NAME}_1.fastq.gz"
    INPUT_READ_FILE_2="${READ_DIR}/${SAMPLE_NAME}_2.fastq.gz"
    GRAPH_FILES_DIR_PATH="/data/markellocj/graph_reference/snp1kg_256-HS37D5"
    WORK_DIR="/data/Udpbinfo/usr/markellocj/${COHORT_NAME}_run/output_${SAMPLE_NAME}_vg_t148"
    VG_CONTAINER="v1.6.0-939-g9fa94410-t149-run"
    mkdir -p ${WORK_DIR} && cd ${WORK_DIR}
    sample_master_script.sh -i ${SAMPLE_NAME} -r ${INPUT_READ_FILE_1} -r ${INPUT_READ_FILE_2} -g ${GRAPH_FILES_DIR_PATH} -w ${WORK_DIR} -c ${VG_CONTAINER} > ${WORK_DIR}/master_script.stdout
    
### Things to keep in mind. Basic job maintenance

Running a pipeline this large can be unstable due to the amount of data and the number of swarm jobs it executes. It's best to keep track of your running jobs by periodically running the following monitoring tools:
- `squeue -u ${USERNAME}` OR `jobload -u ${USERNAME}`: Gives a list of currently running slurm jobs on the system along with JOBID's and current resource use for each job.
- `jobhist ${JOBID}`: Gives specific execution time and success status of previously execulted jobs.
- `checkquota -u ${USERNAME}`: Gives total disk quota provided and disk useage for the user.

More information on monitorying Slurm batchsystem jobs on Biowulf can be found [here](https://hpc.nih.gov/docs/userguide.html#monitor).

### Run times

The mapping pipeline takes approximately between 5-10 hours, but can very significantly based on the availability of nodes. You can get a list of available nodes by using the `freen` command.

