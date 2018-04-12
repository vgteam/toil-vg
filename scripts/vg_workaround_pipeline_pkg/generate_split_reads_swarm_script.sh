#!/bin/bash

# Split reads
# Assumes that source ead files are located in the same directory

# Resource recommendation for running this script
#   sbatch --cpus-per-task=32 --mem=200g --time=30:00:00
# Should take about 30 minutes 

# Example run_split_reads.sh /data/Udpbinfo/usr/markellocj/sandbox/ERR174310_1.fastq.gz /data/Udpbinfo/usr/markellocj/sandbox/ERR174310_2.fastq.gz NA12877 /data/Udpbinfo/usr/markellocj/sandbox quay.io/vgteam/vg:v1.5.0-2018-g71f96239-t119-run

# Example run_split_reads.sh /data/Udpdata/Individuals/UDP10618/WGS/Rawreads/HudsonAlpha_1/UDN021466-H3CGJALXX_s1_1_GSLv3-7_01_SL202061.fastq.gz /data/Udpdata/Individuals/UDP10618/WGS/Rawreads/HudsonAlpha_1/UDN021466-H3CGJALXX_s1_2_GSLv3-7_01_SL202061.fastq.gz UDP10618 /data/Udpbinfo/usr/markellocj/sandbox quay.io/vgteam/vg:v1.5.0-2018-g71f96239-t119-run

INPUT_READ_FILE_1=$1
INPUT_READ_FILE_2=$2
SAMPLE_NAME=$3
WORK_DIR=$4
VG_CONTAINER=$5
FASTQ_SPLIT_SWARMFILE_NAME="${WORK_DIR}/split_fastq_swarmfile_${SAMPLE_NAME}"
FASTQ_WORKDIR_PATH="${WORK_DIR}/split_fastqs_${SAMPLE_NAME}"

rm ${FASTQ_SPLIT_SWARMFILE_NAME}

## Check if the fastq working directory exists for the user running this script; if not, make the directory
if [ ! -d "$FASTQ_WORKDIR_PATH" ]; then
    mkdir -p $FASTQ_WORKDIR_PATH
    chmod 2770 $FASTQ_WORKDIR_PATH
fi

READS_PER_CHUNK=10000000
let CHUNK_LINES=$READS_PER_CHUNK*4

echo "module load singularity; cd ${FASTQ_WORKDIR_PATH}; gzip -cd ${INPUT_READ_FILE_1} | SINGULARITYENV_CHUNK_LINES=${CHUNK_LINES} singularity exec -H ${FASTQ_WORKDIR_PATH}:${HOME} --pwd ${HOME} docker://${VG_CONTAINER} split -l ${CHUNK_LINES} --filter='pigz -p 32 > \$FILE.fq.gz' - ${HOME}/fq_chunk_1.part." >> ${FASTQ_SPLIT_SWARMFILE_NAME}

echo "module load singularity; cd ${FASTQ_WORKDIR_PATH}; gzip -cd ${INPUT_READ_FILE_2} | SINGULARITYENV_CHUNK_LINES=${CHUNK_LINES} singularity exec -H ${FASTQ_WORKDIR_PATH}:${HOME} --pwd ${HOME} docker://${VG_CONTAINER} split -l ${CHUNK_LINES} --filter='pigz -p 32 > \$FILE.fq.gz' - ${HOME}/fq_chunk_2.part." >> ${FASTQ_SPLIT_SWARMFILE_NAME}

