#!/bin/bash


SAMPLE_NAME=$1
WORK_DIR=$2
GRAPH_FILES_DIR_PATH=$3
ID_RANGES_FILE_NAME="genome_id_ranges.tsv"
FASTQ_CHUNK_DIR="${WORK_DIR}/split_fastqs_${SAMPLE_NAME}"
BED_DIR_NAME="gam_refactor_bed_${SAMPLE_NAME}"
BED_DIR_PATH="${WORK_DIR}/${BED_DIR_NAME}"
GAM_DIR_NAME="swarm_map_gam_output_${SAMPLE_NAME}"
GAM_DIR_PATH="${WORK_DIR}/${GAM_DIR_NAME}"
CHR_OUTPUT_GAM_DIR_NAME="${SAMPLE_NAME}_final_output_gams"
CHR_OUTPUT_GAM_DIR_PATH="${WORK_DIR}/${CHR_OUTPUT_GAM_DIR_NAME}"


cd ${WORK_DIR}

if [ ! -d "${CHR_OUTPUT_GAM_DIR_PATH}" ]; then
    mkdir -p ${CHR_OUTPUT_GAM_DIR_PATH}
    chmod 2770 ${CHR_OUTPUT_GAM_DIR_PATH}
fi

read_chunk_ids=($(ls -l $FASTQ_CHUNK_DIR | awk -F'.' '{print $3}' | sort | uniq | xargs))

chr_id_list=()
index_num=0

while IFS='' read -r line; do
    chr_id=($(echo "$line" | cut -d$'\t' -f1))
    chr_id_list+=(${chr_id})
    index_num=$(($index_num + 1))
done < "${GRAPH_FILES_DIR_PATH}/${ID_RANGES_FILE_NAME}"

## Extract filenames from chunked .bed files

# For each chromosome create a list which will each contain the chunks for that particular chromosome

declare -A chunked_chr_id_file_list=()

# For each chunk id
for i in ${read_chunk_ids[@]}
do
    BED_FILE_NAME="${BED_DIR_PATH}/${i}_output_bed.bed"
    index_num=0
    while IFS='' read -r line; do
        chr_id=${chr_id_list[$(($index_num))]}
        chunked_chr_id_file_list[$chr_id]+=" ${GAM_DIR_PATH}/$(echo "$line" | cut -d$'\t' -f4 | rev | cut -d$'/' -f1 | rev)"
        index_num=$(($index_num + 1))
    done < "${BED_FILE_NAME}"
done

# For each chromosome, concatenate the chromosome chunked .gams into a single chromosomal .gam file
for j in ${chr_id_list[@]}
do
    cat ${chunked_chr_id_file_list[$j]} > "${CHR_OUTPUT_GAM_DIR_PATH}/${SAMPLE_NAME}_${j}.gam"
done

