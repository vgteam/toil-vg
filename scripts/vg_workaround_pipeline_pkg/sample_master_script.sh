#!/bin/bash
#################################################################################################
##
##  Script to run the full VG workaround pipeline (up to surjected BAMs) for Whole Genome data
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
if [ $# -lt 5 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

PACKAGE_DIR="/data/markellocj/vg_workaround_pipeline_pkg"

VG_CONTAINER="quay.io/vgteam/vg:v1.6.0-187-ga5bc5549-t124-run"
## Parse through arguments
while getopts "i:r:g:w:c:h" OPTION; do
    case $OPTION in
        i)
            SAMPLE_NAME=$OPTARG
        ;;
        r)
            FASTQ_FILES+=($OPTARG)
        ;;
        g)
            GRAPH_FILES_DIR=$OPTARG
        ;;
        w)
            WORK_DIR=$OPTARG
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

## STEP1: SPLIT READS. Generate and run read splitting swarm jobscript.
INPUT_READ_FILE_1=${FASTQ_FILES[0]}
INPUT_READ_FILE_2=${FASTQ_FILES[1]}

echo "SAMPLE_NAME: ${SAMPLE_NAME}"
echo "INPUT_READ_FILE_1: ${INPUT_READ_FILE_1}"
echo "INPUT_READ_FILE_2: ${INPUT_READ_FILE_2}"
echo "GRAPH_FILES_DIR: ${GRAPH_FILES_DIR}"
echo "WORK_DIR: ${WORK_DIR}"
echo "VG_CONTAINER: ${VG_CONTAINER}"

if [ ! -d "${WORK_DIR}" ]; then
    mkdir -p ${WORK_DIR}
    chmod 2770 ${WORK_DIR}
fi

echo "Generating swarm script for splitting reads."
${PACKAGE_DIR}/generate_split_reads_swarm_script.sh ${INPUT_READ_FILE_1} ${INPUT_READ_FILE_2} ${SAMPLE_NAME} ${WORK_DIR} ${VG_CONTAINER}
echo "Finished generating swarm script for splitting reads."

SPLIT_READS_SWARMFILE_PATH="${WORK_DIR}/split_fastq_swarmfile_${SAMPLE_NAME}"
SPLIT_READS_JOBID=$(swarm -f ${SPLIT_READS_SWARMFILE_PATH} -g 100 -t 32 --time=2:00:00)
echo "Splitting reads into chunks. Jobid: ${SPLIT_READS_JOBID}"

## STEP2: GENERATE SWARM SCRIPTS. Generate the swarm scripts used for running the chunk alignment,
##        chromosome split, and indexing chromosomal gams procedures.
GENERATE_SWARM_JOBID=$(sbatch --time=1:00:00 --dependency=afterok:${SPLIT_READS_JOBID} ${PACKAGE_DIR}/generate_vg_map_pipeline_swarm_scripts.sh ${SAMPLE_NAME} ${WORK_DIR} ${GRAPH_FILES_DIR} ${VG_CONTAINER} ${PACKAGE_DIR})
echo "Generating swarm scripts for alignment pipeline. Jobid:${GENERATE_SWARM_JOBID}"

## STEPS 3-7 (CHUNK ALIGNMENT, INDEX GAMs, and CHROMOSOME SPLIT, MERGE GAMs, SURJECT) are run within STEP2 to workaround
##        pipeline dependencies for dynamically generated swarm scripts.

exit

