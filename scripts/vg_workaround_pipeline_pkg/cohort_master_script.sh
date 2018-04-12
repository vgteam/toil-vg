#!/bin/bash
#################################################################################################
##
##  Script to run an entire UDP cohort through the full VG workaround pipeline
##      (up to surjected BAMs) for Whole Genome data.
##
##  Inputs:
##
##  Assumptions:
##      - Reads are located in a different directory than in the working directory of where
##          the argument WORK_DIR is set.
##          - This is to get around the problem of binding two different directories when
##              running Singularity containers.
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script runs the VG pipeline up to GAM surjection to BAM files (no variant calling is made).

Inputs:
    -i Sample ID (in format UDP####)

Outputs:

Assumptions:

EOF

}

## Check whether script is being run on Biowulf
if [ "$HOSTNAME" == helix.nih.gov ]; then
    usage
    exit 1
fi

## Check number of arguments
if [ $# -lt 4 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

PACKAGE_DIR="/data/Udpbinfo/usr/markellocj/vg_workaround_pipeline_pkg"

VG_CONTAINER="quay.io/vgteam/vg:v1.6.0-187-ga5bc5549-t124-run"
## Parse through arguments
while getopts "i:g:w:c:h" OPTION; do
    case $OPTION in
        i)
            COHORT_NAME=$OPTARG
        ;;
        g)
            GRAPH_FILES_DIR=$OPTARG
        ;;
        w)
            COHORT_SET_WORK_DIR=$OPTARG
        ;;
        c)
            VG_CONTAINER=$OPTARG
        ;;
        h)
            usage
            exit 1
        ;;
        \?)
            usage
            exit 1
        ;;
    esac
done

## STEP1: COLLECT READS. Generate sample work directories and create read gathering swarm script
COHORT_WORK_DIR="${COHORT_SET_WORK_DIR}/${COHORT_NAME}_run"
COHORT_DATA_DIR="/data/Udpdata/Families/${COHORT_NAME}"
COHORT_NAMES_LIST=($(ls $COHORT_DATA_DIR | grep 'UDP'))
CAT_READS_SWARMFILE_PATH="${COHORT_WORK_DIR}/cat_reads_swarmfile_${COHORT_NAME}"

mkdir ${COHORT_WORK_DIR}

echo "Cohort UDP Sample List: ${COHORT_NAMES_LIST[@]}"
touch cat_reads_swarm_script

for SAMPLE_NAME in ${COHORT_NAMES_LIST[@]}
do
    ## Make readfile directories
    mkdir ${COHORT_WORK_DIR}/${SAMPLE_NAME}_reads

    ## Concatenate reads and put them in the sample-specific readfile directory
    PAIR_1_READS=()
    PAIR_2_READS=()
    LANE_NUMS=($(ls ${COHORT_DATA_DIR}/${SAMPLE_NAME}/WGS/Rawreads/HudsonAlpha_1 | awk -F'_' '{print $2}' | sort | uniq | xargs))
    for LANE_NUM in ${LANE_NUMS[@]}
    do
        PAIR_1_READS+=(${COHORT_DATA_DIR}/${SAMPLE_NAME}/WGS/Rawreads/HudsonAlpha_1/"$(ls ${COHORT_DATA_DIR}/${SAMPLE_NAME}/WGS/Rawreads/HudsonAlpha_1 | grep "_${LANE_NUM}_1")")
        PAIR_2_READS+=(${COHORT_DATA_DIR}/${SAMPLE_NAME}/WGS/Rawreads/HudsonAlpha_1/"$(ls ${COHORT_DATA_DIR}/${SAMPLE_NAME}/WGS/Rawreads/HudsonAlpha_1 | grep "_${LANE_NUM}_2")")
    done

    echo "cat ${PAIR_1_READS[@]} > ${COHORT_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_1.fq.gz" >> ${CAT_READS_SWARMFILE_PATH}
    echo "cat ${PAIR_2_READS[@]} > ${COHORT_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_2.fq.gz" >> ${CAT_READS_SWARMFILE_PATH}
done

COLLECT_READS_JOBID=$(swarm -f ${CAT_READS_SWARMFILE_PATH} --time=4:00:00)
echo "Running read collection swarm script. Jobid:${COLLECT_READS_JOBID}"

cd ${COHORT_WORK_DIR}

## STEP2: RUN MAPPING PIPELINE. Run the graph mapping pipeline per sample.
for SAMPLE_NAME in ${COHORT_NAMES_LIST[@]}
do
    READ_FILE_1="${COHORT_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_1.fq.gz"
    READ_FILE_2="${COHORT_WORK_DIR}/${SAMPLE_NAME}_reads/${SAMPLE_NAME}_read_pair_2.fq.gz"
    SAMPLE_MAP_RUN_JOBID=$(sbatch --time=12:00:00 --dependency=afterok:${COLLECT_READS_JOBID} ${PACKAGE_DIR}/sample_master_script.sh -i ${SAMPLE_NAME} -r ${READ_FILE_1} -r ${READ_FILE_2} -g ${GRAPH_FILES_DIR} -w ${COHORT_WORK_DIR} -c "quay.io/vgteam/vg:v1.6.0-187-ga5bc5549-t124-run")
    echo "Running sample ${SAMPLE_NAME} through the graph alignment pipeline. Jobid:${SAMPLE_MAP_RUN_JOBID}"
done


exit

