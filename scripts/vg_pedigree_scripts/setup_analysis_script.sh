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
    -f Cohort name. Should be the same as the proband sample name.
    -c PATH to .ped file containing all samples in the family.
    -w PATH to where the UDP cohort will be processed.
    -a PATH to chromosome annotation directory used by vcftoshebang.
    -e PATH to directory containing master edit files used by vcftoshebang.
    -d PATH to cadd engine data directory.
    -v PATH to the toil_vg repository
    -r (OPTIONAL, default=false) Set to 'true' to restart an incompletely ran workflow
    
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
RESTART=false

## Parse through arguments
while getopts "f:c:w:a:e:d:v:r:h" OPTION; do
    case $OPTION in
        f)
            COHORT_NAME=$OPTARG
        ;;
        c)
            COHORT_PED_FILE=$OPTARG
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
        ;;
        a)
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

# Extract sample information from input family .ped file
source ${TOIL_VG_DIR}/toilvg_venv/bin/activate
pip install ped_parser
TRIO_PED_FILE="${COHORT_NAME}.trio.ped"

READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
SAMPLES_LIST=($(python3 -c "import ped_parser; print(list(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals.keys()))" | tr -d '[],' | tr -d \'))
SIB_ID_LIST_SET=()
SIB_GENDER_LIST_SET=()
SIB_AFFECTED_LIST_SET=()

SIBLING_SAMPLE_NAMES=()
SIBLING_GENDERS=()
SIBLING_AFFECTED=()

for SAMPLE_ID in ${SAMPLES_LIST[@]}
do
    SAMPLE_MOM=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].mother)"))
    SAMPLE_DAD=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].father)"))
    SAMPLE_GENDER=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].sex)"))
    SAMPLE_AFFECTED=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].phenotype)"))
    if [[ "$SAMPLE_ID" == "$COHORT_NAME" ]]; then
        SIBLING_SAMPLE_NAMES+=(${SAMPLE_ID})
        SIBLING_GENDERS+=($((${SAMPLE_GENDER} - 1)))
        if [[ ${SAMPLE_AFFECTED} -eq 1 ]]; then
            SIBLING_AFFECTED+=(1)
        else
            SIBLING_AFFECTED+=(0)
        fi
        MATERNAL_SAMPLE_NAME="${SAMPLE_MOM}"
        PATERNAL_SAMPLE_NAME="${SAMPLE_DAD}"
    elif [[ ${SAMPLE_MOM} != '0' ]]; then
        SIB_ID_LIST_SET+=(${SAMPLE_ID})
        SIB_GENDER_LIST_SET+=(${SAMPLE_GENDER})
        SIB_AFFECTED_LIST_SET+=(${SAMPLE_AFFECTED})
    fi
done

for (( n=0; n<${#SIB_ID_LIST_SET[@]}; n++ ))
do
    SIBLING_SAMPLE_NAMES+=(${SIB_ID_LIST_SET[$n]})
    SIBLING_GENDERS+=($((${SIB_GENDER_LIST_SET[$n]} - 1)))
    if [[ ${SIB_AFFECTED_LIST_SET[$n]} -eq 1 ]]; then
        SIBLING_AFFECTED+=(1)
    else
        SIBLING_AFFECTED+=(0)
    fi
done

# Format the input values in the command by the parsed sample information
PROBAND_SAMPLE_NAME="${SIBLING_SAMPLE_NAMES[0]}"
INPUT_DATA_DIR="${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore"
SIB_ID_LIST=""
SIB_GENDER_LIST=""
SIB_AFFECTED_LIST=""
SIB_BAM_LIST=""
SIB_BAI_LIST=""
for (( n=0; n<${#SIBLING_SAMPLE_NAMES[@]}; n++ ))
do
    SIB_ID_LIST+="${SIBLING_SAMPLE_NAMES[$n]} "
    SIB_GENDER_LIST+="${SIBLING_GENDERS[$n]} "
    SIB_AFFECTED_LIST+="${SIBLING_AFFECTED[$n]} "
    SIB_BAM_LIST+="'${INPUT_DATA_DIR}/${SIBLING_SAMPLE_NAMES[$n]}_merged.indel_realigned.bam' "
    SIB_BAI_LIST+="'${INPUT_DATA_DIR}/${SIBLING_SAMPLE_NAMES[$n]}_merged.indel_realigned.bai' "
done

if [[ ${COHORT_WORKFLOW_DIR} = *[[:space:]]* ]] || [ -z ${COHORT_WORKFLOW_DIR} ]; then
    echo "ERROR: ${COHORT_WORKFLOW_DIR} argument value contains whitespace or is empty"
    exit 1
fi
if [[ ${PROBAND_SAMPLE_NAME} = *[[:space:]]* ]] || [ -z ${PROBAND_SAMPLE_NAME} ]; then
    echo "ERROR: ${PROBAND_SAMPLE_NAME} argument value contains whitespace or is empty"
    exit 1
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
--statePollingWait 60 \\
--rescueJobsFrequency 60 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir onSuccess \\
--whole_genome_config \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_jobstore \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore \\
--cohort_vcf ${INPUT_DATA_DIR}/${PROBAND_SAMPLE_NAME}.snpeff.unrolled.vcf.gz \\
--sample_name ${PROBAND_SAMPLE_NAME} \\
--maternal_name ${MATERNAL_SAMPLE_NAME} \\
--paternal_name ${PATERNAL_SAMPLE_NAME} \\
--sibling_names ${SIB_ID_LIST[@]} \\
--sibling_genders ${SIB_GENDER_LIST[@]} \\
--sibling_affected ${SIB_AFFECTED_LIST[@]} \\
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

