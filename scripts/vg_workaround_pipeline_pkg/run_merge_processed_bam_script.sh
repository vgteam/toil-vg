#!/bin/bash

SAMPLE_NAME=$1
WORK_DIR=$2
CHR_OUTPUT_BAM_DIR_NAME="${SAMPLE_NAME}_surjected_chr_bams"
CHR_OUTPUT_BAM_DIR_PATH="${WORK_DIR}/${CHR_OUTPUT_BAM_DIR_NAME}"
PROCESS_BAM_CORES=32

cd ${WORK_DIR} 

# Merge the chromosomal BAM files into a single sample wgs BAM file
module load samtools
samtools merge -f --threads ${PROCESS_BAM_CORES} ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_merged.dupmarked.reordered.bam ${CHR_OUTPUT_BAM_DIR_PATH}/*.sorted.dupmarked.reordered.bam

# Index the final merged and reordered BAM file
samtools index -@ ${PROCESS_BAM_CORES} ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_merged.dupmarked.reordered.bam

