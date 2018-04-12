#!/bin/bash

## Assumes xg, gcsa, id_ranges, and cid_ranges files names

SAMPLE_NAME=$1
WORK_DIR=$2
GRAPH_FILES_DIR_PATH=$3
VG_CONTAINER=$4
PACKAGE_DIR=$5
XG_FILE_NAME="1kg_ref_hs37d5.xg"
GCSA_FILE_NAME="1kg_ref_hs37d5.gcsa"
ID_RANGES_FILE_NAME="genome_id_ranges.tsv"
CID_RANGES_FILE_NAME="cid_ranges.txt"
FASTQ_DIR_NAME="split_fastqs_${SAMPLE_NAME}"
SWARM_MAP_OUTPUT_DIR_NAME="swarm_map_gam_output_${SAMPLE_NAME}"
BED_DIR_NAME="gam_refactor_bed_${SAMPLE_NAME}"

BED_WORKDIR_PATH="${WORK_DIR}/${BED_DIR_NAME}"
BED_CONTAINER_PATH="${HOME}/${BED_DIR_NAME}"
SWARM_MAP_OUTPUT_WORKDIR_PATH="${WORK_DIR}/${SWARM_MAP_OUTPUT_DIR_NAME}"
SWARM_MAP_OUTPUT_CONTAINER_PATH="${HOME}/${SWARM_MAP_OUTPUT_DIR_NAME}"
FASTQ_CHUNK_WORKDIR_PATH="${WORK_DIR}/${FASTQ_DIR_NAME}"

MAP_SWARMFILE_NAME="${WORK_DIR}/map_swarmfile_${SAMPLE_NAME}"
INDEX_SWARMFILE_NAME="${WORK_DIR}/index_gams_chrom_swarmfile_${SAMPLE_NAME}"
SPLIT_GAMS_SWARMFILE_NAME="${WORK_DIR}/split_gams_chrom_swarmfile_${SAMPLE_NAME}"
ALIGNMENT_CORES=32
INDEX_CORES=32
SPLIT_CORES=32

# Extract chunk ids (e.g. fq_chunk_1.part.aa.fq.gz and fq_chunk_2.part.aa.fq.gz --> aa)
read_chunk_ids=($(ls -l ${FASTQ_CHUNK_WORKDIR_PATH} | awk -F'.' '{print $3}' | sort | uniq | xargs))

rm ${MAP_SWARMFILE_NAME} ${INDEX_SWARMFILE_NAME} ${SPLIT_GAMS_SWARMFILE_NAME} ${GRAPH_FILES_DIR_PATH}/${CID_RANGES_FILE_NAME}

if [ ! -d "${SWARM_MAP_OUTPUT_WORKDIR_PATH}" ]; then
    mkdir -p ${SWARM_MAP_OUTPUT_WORKDIR_PATH}
    chmod 2770 ${SWARM_MAP_OUTPUT_WORKDIR_PATH}
fi

if [ ! -d "${BED_WORKDIR_PATH}" ]; then
    mkdir -p ${BED_WORKDIR_PATH}
    chmod 2770 ${BED_WORKDIR_PATH}
fi

while IFS='' read -r line; do
    echo "$line" | cut -d$'\t' -f2,3 --output-delimiter=":" >> "${GRAPH_FILES_DIR_PATH}/${CID_RANGES_FILE_NAME}"
done < "${GRAPH_FILES_DIR_PATH}/${ID_RANGES_FILE_NAME}"


# Iterate for each chunk id and run vg map for each chunked fastq
for i in ${read_chunk_ids[@]}
do
    echo $i
    READ_FILE_1_CONTAINER_PATH=${HOME}/${FASTQ_DIR_NAME}/fq_chunk_1.part.${i}.fq.gz
    READ_FILE_2_CONTAINER_PATH=${HOME}/${FASTQ_DIR_NAME}/fq_chunk_2.part.${i}.fq.gz
    map_command="module load singularity; cd ${WORK_DIR}; singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} -B ${GRAPH_FILES_DIR_PATH}:/mnt docker://${VG_CONTAINER} vg mpmap -S -f ${READ_FILE_1_CONTAINER_PATH} -f ${READ_FILE_2_CONTAINER_PATH} -x /mnt/${XG_FILE_NAME} -g /mnt/${GCSA_FILE_NAME} -t ${ALIGNMENT_CORES} > ${SWARM_MAP_OUTPUT_WORKDIR_PATH}/${SAMPLE_NAME}_${i}.gam"
    #map_command="module load singularity/2.4; cd ${WORK_DIR}; singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} -B ${GRAPH_FILES_DIR_PATH}:/mnt docker://${VG_CONTAINER} vg map --read-group "${SAMPLE_NAME}_reads" --sample ${SAMPLE_NAME} -f ${READ_FILE_1_CONTAINER_PATH} -f ${READ_FILE_2_CONTAINER_PATH} -x /mnt/${XG_FILE_NAME} -g /mnt/${GCSA_FILE_NAME} -t ${ALIGNMENT_CORES} > ${SWARM_MAP_OUTPUT_WORKDIR_PATH}/${SAMPLE_NAME}_${i}.gam"

    echo $map_command >> ${MAP_SWARMFILE_NAME}
    
    OUTPUT_BED_NAME="${BED_WORKDIR_PATH}/${i}_output_bed.bed"
    touch ${OUTPUT_BED_NAME}

    index_command="module load singularity; cd ${WORK_DIR}; singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} docker://${VG_CONTAINER} vg index -a ${SWARM_MAP_OUTPUT_CONTAINER_PATH}/${SAMPLE_NAME}_${i}.gam -d ${SWARM_MAP_OUTPUT_CONTAINER_PATH}/${SAMPLE_NAME}_${i}.gam.index -t ${INDEX_CORES}"

    echo $index_command >> ${INDEX_SWARMFILE_NAME}

    split_command="module load singularity; cd ${WORK_DIR}; singularity -q exec -H ${WORK_DIR}:${HOME} --pwd ${HOME} -B ${GRAPH_FILES_DIR_PATH}:/mnt docker://${VG_CONTAINER} vg chunk -x /mnt/${XG_FILE_NAME} -a ${SWARM_MAP_OUTPUT_CONTAINER_PATH}/${SAMPLE_NAME}_${i}.gam.index -c 0 -R /mnt/${CID_RANGES_FILE_NAME} -b ${SWARM_MAP_OUTPUT_CONTAINER_PATH}/${SAMPLE_NAME}_${i}.gam -t ${SPLIT_CORES} -E ${BED_CONTAINER_PATH}/${i}_output_bed.bed"

    echo $split_command >> ${SPLIT_GAMS_SWARMFILE_NAME}
done


## STEP3: CHUNK ALIGNMENT.
CHUNK_ALIGNMENT_JOBID=$(swarm -f ${MAP_SWARMFILE_NAME} -g 100 -t 32 --time 8:00:00 --maxrunning 15)
echo "Running swarm VG graph alignment. Jobid:${CHUNK_ALIGNMENT_JOBID}"

## STEP4: INDEX GAMs.
INDEX_GAM_JOBID=$(swarm -f ${INDEX_SWARMFILE_NAME} -g 100 -t 32 --time 1:00:00 --maxrunning 15 --dependency=afterok:${CHUNK_ALIGNMENT_JOBID})
echo "Running swarm gam indexing. Jobid:${INDEX_GAM_JOBID}"

## STEP5: CHROMOSOME SPLIT.
SPLIT_GAMS_JOBID=$(swarm -f ${SPLIT_GAMS_SWARMFILE_NAME} -g 100 -t 32 --time 2:00:00 --maxrunning 15 --dependency=afterok:${INDEX_GAM_JOBID})
echo "Running swarm script to split GAMs by chromosome. Jobid:${SPLIT_GAMS_JOBID}"

## STEP6: MERGE GAMS.
MERGE_GAMS_JOBID=$(sbatch --cpus-per-task=32 --mem=100g --time=1:00:00 --dependency=afterok:${SPLIT_GAMS_JOBID} ${PACKAGE_DIR}/run_merge_chrom_gam_script.sh ${SAMPLE_NAME} ${WORK_DIR} ${GRAPH_FILES_DIR_PATH})
echo "Running merge GAMs into chromosomal GAMs. Jobid:${MERGE_GAMS_JOBID}"

## STEP7: SURJECT. Generate swarm script to run vg surject on each chromosome GAM into final BAM files.
GENERATE_SURJECT_SWARM_JOBID=$(sbatch --time=1:00:00 --dependency=afterok:${MERGE_GAMS_JOBID} ${PACKAGE_DIR}/generate_vg_surject_swarm_script.sh ${SAMPLE_NAME} ${WORK_DIR} ${GRAPH_FILES_DIR_PATH} ${VG_CONTAINER} ${PACKAGE_DIR})
echo "Generating swarm script for surjecting GAM files into BAM files. Jobid:${GENERATE_SURJECT_SWARM_JOBID}"


