#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the entire toil-vg candidate analysis pipeline.
##
##  Inputs:
##
##  Assumptions:
##      The UDP cohort ID is the same as the UDP ID for the proband sample.
##      The working directory as specified by the -w argument is the same as the working
##          directory for the run of 'toil-vg pedigree'
##  
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script setups up a bash script to run a UDP cohort through all stages of the toil-vg
candidate analysis workflow on the NIH Biowulf Cluster.

Inputs:
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -s List of Sibling UDP ID, Proband ID must be first in the list (in format UDP#### UDP#### UDP#### ...)
    -g List of Sibling Gender IDs. 0=male, 1=female. Must be same order as the input to -s argument.
    -a List of Sibling affected status. 0=unaffected, 1=affected. Must be same order as the input to -s argument.
    -w PATH to where the UDP cohort will be processed.
    -c PATH to chromosome annotation directory used by vcftoshebang.
    -e PATH to directory containing master edit files used by vcftoshebang.
    -d PATH to cadd engine data directory.
    -v PATH to the toil_vg repository
    -r (OPTIONAL, default=false) Set to 'true' to restart an incompletely ran workflow
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 7 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
RUN_SMALL_TEST=false
RESTART=false

## Parse through arguments
while getopts "m:f:s:w:i:c:e:d:v:r:t:h" OPTION; do
    case $OPTION in
        m)
            MATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        f)
            PATERNAL_SAMPLE_NAME=$OPTARG
        ;;
        s)
            SIBLING_SAMPLE_NAMES+=($OPTARG)
        ;;
        g)
            SIBLING_GENDERS+=($OPTARG)
        ;;
        a)
            SIBLING_AFFECTED+=($OPTARG)
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
        ;;
        c)
            CHROM_ANNOT_DIR=$OPTARG
        ;;
        e)
            EDIT_ANNOT_DIR=$OPTARG
        ;;
        d)
            CADD_DATA_DIR=$OPTARG
        ;;
        v)
            TOIL_VG_DIR=$OPTARG
        ;;
        r)
            RESTART=$OPTARG
        ;;
        t)
            RUN_SMALL_TEST=$OPTARG
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

PROBAND_SAMPLE_NAME="${SIBLING_SAMPLE_NAMES[0]}"
INPUT_DATA_DIR="${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore"
SIB_ID_LIST=""
SIB_GENDER_LIST=""
SIB_AFFECT_LIST=""
SIB_BAM_LIST=""
SIB_BAI_LIST=""
for (( n=1; n<=${#SIBLING_SAMPLE_NAMES[@]}; n++ ))
do
    SIB_ID_LIST+="'${SIBLING_SAMPLE_NAMES[$n]}' "
    SIB_GENDER_LIST+="'${SIBLING_GENDERS[$n]}' "
    SIB_AFFECT_LIST+="'${SIBLING_AFFECTED[$n]}' "
    SIB_BAM_LIST+="'${INPUT_DATA_DIR}/${SIBLING_SAMPLE_NAMES[$n]}_merged.indel_realigned.bam' "
    SIB_BAI_LIST+="'${INPUT_DATA_DIR}/${SIBLING_SAMPLE_NAMES[$n]}_merged.indel_realigned.bai' "
done

if [[ ${COHORT_WORKFLOW_DIR} = *[[:space:]]* ]]; then
    echo "ERROR: ${COHORT_WORKFLOW_DIR} argument value contains whitespace"
    exit 1
fi
if [[ ${PROBAND_SAMPLE_NAME} = *[[:space:]]* ]]; then
    echo "ERROR: ${PROBAND_SAMPLE_NAME} argument value contains whitespace"
    exit 1
fi

if [ $RESTART == false ]; then
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore
    rm -fr ${COHORT_WORKFLOW_DIR}/tmp
fi

if [ ! -d "${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore"
fi
if [ ! -d "${COHORT_WORKFLOW_DIR}/tmp" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/tmp"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/tmp"
fi

if [ ! -d "/data/$USER/singularity_cache" ]; then
    mkdir -p "/data/$USER/singularity_cache"
    chmod 2770 "/data/$USER/singularity_cache"
fi

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
echo "module load singularity python/3.7" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
echo "source ${TOIL_VG_DIR}/toilvg_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
echo "export TOIL_SLURM_ARGS='-t 20:00:00'" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
echo "export SINGULARITY_CACHEDIR=/data/$USER/singularity_cache" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
if [ $RESTART == false ]; then
    echo "toil clean ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_jobstore" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
fi
RESTART_ARG=""
if [ $RESTART == true ]; then
    RESTART_ARG="--restart"
fi
echo "toil-vg analysis \\
${RESTART_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 30 \\
--rescueJobsFrequency 30 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir onSuccess \\
--whole_genome_config \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_jobstore \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore \\
--cohort_vcf ${INPUT_DATA_DIR}/${PROBAND_SAMPLE_NAME}.snpeff.unrolled.vcf \\
--sample_name ${PROBAND_SAMPLE_NAME} \\
--maternal_name ${MATERNAL_SAMPLE_NAME} \\
--paternal_name ${PATERNAL_SAMPLE_NAME} \\
--sibling_names ${SIB_ID_LIST[@]} \\
--sibling_genders ${SIB_GENDER_LIST[@]} \\
--sibling_affected ${SIB_AFFECT_LIST[@]} \\
--maternal_bam ${INPUT_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam \\
--maternal_bai ${INPUT_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam.bai \\
--paternal_bam ${INPUT_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam \\
--paternal_bai ${INPUT_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_merged.indel_realigned.bam.bai \\
--siblings_bam ${SIB_BAM_LIST[@]} \\
--siblings_bai ${SIB_BAI_LIST[@]} \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh

exit

