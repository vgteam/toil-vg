#!/bin/bash

SAMPLE_NAME=$1
WORK_DIR=$2
PACKAGE_DIR=$3
SPLIT_FASTQ_DIR_NAME="split_fastqs_${SAMPLE_NAME}"
GAM_REFACTOR_BED_DIR_NAME="gam_refactor_bed_${SAMPLE_NAME}"
MAP_GAM_SWARM_OUTPUT_DIR_NAME="swarm_map_gam_output_${SAMPLE_NAME}"

SPLIT_FASTQ_SWARMFILE_NAME="split_fastq_swarmfile_${SAMPLE_NAME}"
INDEX_GAMS_SWARMFILE_NAME="index_gams_chrom_swarmfile_${SAMPLE_NAME}"
MAP_SWARMFILE_NAME="map_swarmfile_${SAMPLE_NAME}"
SPLIT_GAMS_CHROM_SWARMFILE_NAME="split_gams_chrom_swarmfile_${SAMPLE_NAME}"
SURJECT_GAMS_SWARMFILE_NAME="surject_gams_swarmfile_${SAMPLE_NAME}"
PROCESS_BAMS_SWARMFILE_NAME="process_bams_swarmfile_${SAMPLE_NAME}"


cd ${WORK_DIR}

# Clean out temporary data directories
rm -fr ${GAM_REFACTOR_BED_DIR_NAME} ${MAP_GAM_SWARM_OUTPUT_DIR_NAME}

# Clean out leftover swarm files
rm -f ${SPLIT_FASTQ_SWARMFILE_NAME} ${INDEX_GAMS_SWARMFILE_NAME} ${MAP_SWARMFILE_NAME} ${SPLIT_GAMS_CHROM_SWARMFILE_NAME} ${SURJECT_GAMS_SWARMFILE_NAME} ${PROCESS_BAMS_SWARMFILE_NAME}

# Clean out swarm stderr and stdout files
cd ${WORK_DIR}; ls | grep '\.o$\|\.e$' | xargs rm

# Clean out temporary surjected bam files
cd "${SAMPLE_NAME}_surjected_chr_bams"; ls | grep '[0-9|Y|X].\(sorted.\|sorted.dupmarked.\|sorted.dupmarked.reordered.\)bam' | xargs rm

