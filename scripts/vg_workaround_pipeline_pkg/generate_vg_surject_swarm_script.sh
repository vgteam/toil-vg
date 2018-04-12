#!/bin/bash

SAMPLE_NAME=$1
WORK_DIR=$2
GRAPH_FILES_DIR_PATH=$3
VG_CONTAINER=$4
PACKAGE_DIR=$5
GAM_DIR_NAME="swarm_map_gam_output_${SAMPLE_NAME}"
GAM_DIR_PATH="${WORK_DIR}/${GAM_DIR_NAME}"
CHR_OUTPUT_GAM_DIR_NAME="${SAMPLE_NAME}_final_output_gams"
CHR_OUTPUT_BAM_DIR_NAME="${SAMPLE_NAME}_surjected_chr_bams"
CHR_OUTPUT_BAM_DIR_PATH="${WORK_DIR}/${CHR_OUTPUT_BAM_DIR_NAME}"
XG_FILE_NAME="1kg_ref_hs37d5.xg"
ID_RANGES_FILE_NAME="genome_id_ranges.tsv"
SURJECT_GAMS_SWARMFILE_NAME="surject_gams_swarmfile_${SAMPLE_NAME}"
SURJECT_CORES=32

cd ${WORK_DIR}

rm -f ${SURJECT_GAMS_SWARMFILE_NAME}

if [ ! -d "${CHR_OUTPUT_BAM_DIR_PATH}" ]; then
    mkdir -p ${CHR_OUTPUT_BAM_DIR_PATH}
    chmod 2770 ${CHR_OUTPUT_BAM_DIR_PATH}
fi

# Clean out old tempoarary chunked GAM directory
rm -fr ${GAM_REFACTOR_BED_DIR_NAME} ${MAP_GAM_SWARM_OUTPUT_DIR_NAME}

# Extract chromosomal path names
CHROMS=""
while IFS='' read -r line; do
    chr_id=($(echo "$line" | cut -d$'\t' -f1))
    CHROMS+=" ${chr_id}"
done < "${GRAPH_FILES_DIR_PATH}/${ID_RANGES_FILE_NAME}"

# Run through each chromosome path name and set a vg surject swarm command
for i in ${CHROMS}
do
    echo $i
    GAM_FILE="${HOME}/${CHR_OUTPUT_GAM_DIR_NAME}/${SAMPLE_NAME}_${i}.gam"
    surject_command="module load singularity; cd ${WORK_DIR}; singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} -B ${GRAPH_FILES_DIR_PATH}:/mnt docker://${VG_CONTAINER} vg surject -x /mnt/${XG_FILE_NAME} -t ${SURJECT_CORES} -b -p ${i} ${GAM_FILE} > ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.bam"
    echo $surject_command >> ${SURJECT_GAMS_SWARMFILE_NAME}
done

# STEP8: RUN SURJECT SWARM SCRIPT
SURJECT_GAMS_JOBID=$(swarm -f ${SURJECT_GAMS_SWARMFILE_NAME} -g 100 -t 32 --time 12:00:00)
echo "Running surject swarm script. Jobid:${SURJECT_GAMS_JOBID}"

# STEP9: SORT, MERGE and REORDER BAMs.
GENERATE_PROCESS_BAMS_SWARM_JOBID=$(sbatch --time=1:00:00 --dependency=afterok:${SURJECT_GAMS_JOBID} ${PACKAGE_DIR}/generate_process_bam_swarm_script.sh ${SAMPLE_NAME} ${WORK_DIR} ${GRAPH_FILES_DIR_PATH} ${PACKAGE_DIR})
echo "Generating swarm script for BAM processing for genome wide variant calling. Jobid:${GENERATE_PROCESS_BAMS_SWARM_JOBID}"

