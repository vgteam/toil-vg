#!/bin/bash

SAMPLE_NAME=$1
WORK_DIR=$2
GRAPH_FILES_DIR_PATH=$3
PACKAGE_DIR=$4
ID_RANGES_FILE_NAME="genome_id_ranges.tsv"
CHR_OUTPUT_BAM_DIR_NAME="${SAMPLE_NAME}_surjected_chr_bams"
CHR_OUTPUT_BAM_DIR_PATH="${WORK_DIR}/${CHR_OUTPUT_BAM_DIR_NAME}"
PROCESS_BAMS_SWARMFILE_NAME="process_bams_swarmfile_${SAMPLE_NAME}"
PROCESS_BAM_CORES=32

cd ${WORK_DIR} 

rm -f ${PROCESS_BAMS_SWARMFILE_NAME}

# Extract chromosomal path names
CHROMS=""
while IFS='' read -r line; do
    chr_id=($(echo "$line" | cut -d$'\t' -f1))
    CHROMS+=" ${chr_id}"
done < "${GRAPH_FILES_DIR_PATH}/${ID_RANGES_FILE_NAME}"

# Run through each chromosome BAM and sort them
# Reorder the merged BAM file to fit the chromosome region order of the reference
for i in ${CHROMS}
do
    echo $i
    process_bam_command="module load samtools picard && samtools sort -T ${CHR_OUTPUT_BAM_DIR_PATH}/tmpSort_${i} --threads ${PROCESS_BAM_CORES} ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.bam > ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.bam && java -Xmx4g -XX:ParallelGCThreads=5 -jar \$PICARDJARPATH/picard.jar MarkDuplicates I=${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.bam O=${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.bam M=${CHR_OUTPUT_BAM_DIR_PATH}/marked_dup_metrics_chr${i}.txt && rm -f ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.bam && java -Xmx20g -XX:ParallelGCThreads=32 -jar \$PICARDJARPATH/picard.jar ReorderSam INPUT=${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.bam OUTPUT=${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.reordered.bam REFERENCE=/data/markellocj/fasta_references/hs37d5_reference/hs37d5.fa && rm -f ${CHR_OUTPUT_BAM_DIR_PATH}/${SAMPLE_NAME}_${i}.sorted.dupmarked.bam"
    echo $process_bam_command >> ${PROCESS_BAMS_SWARMFILE_NAME}
done

PROCESS_BAMS_JOBID=$(swarm -f ${PROCESS_BAMS_SWARMFILE_NAME} -g 100 -t 32 --time 12:00:00 --maxrunning 20)
echo "Running process BAMs swarm script. Jobid:${PROCESS_BAMS_JOBID}"

MERGE_BAMS_JOBID=$(sbatch --cpus-per-task=32 --mem=100g --time=12:00:00 --dependency=afterok:${PROCESS_BAMS_JOBID} ${PACKAGE_DIR}/run_merge_processed_bam_script.sh ${SAMPLE_NAME} ${WORK_DIR})
echo "Running merge of final BAMs. Jobid:${MERGE_BAMS_JOBID}"

# Run the work directory cleanup script if BAM processing completes error free
CLEAN_WORKDIR_JOBID=$(sbatch --time=1:00:00 --dependency=afterok:${MERGE_BAMS_JOBID} ${PACKAGE_DIR}/run_workdir_cleanup.sh ${SAMPLE_NAME} ${WORK_DIR} ${PACKAGE_DIR})
echo "Running work directory cleanup script. Jobid:${CLEAN_WORKDIR_JOBID}"

