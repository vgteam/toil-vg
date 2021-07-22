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
    -b (OPTIONAL, default=true) Set to 'false' to use the GRCh37 cadd annotations
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
GRCh38_REFERENCE_VERSION=true

## Parse through arguments
while getopts "f:c:w:a:e:d:v:r:b:h" OPTION; do
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
        b)
            GRCh38_REFERENCE_VERSION=$OPTARG
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

READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
SAMPLES_LIST=($(python3 -c "import ped_parser; print(list(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals.keys()))" | tr -d '[],' | tr -d \'))
SIB_ID_LIST_SET=()
SIB_GENDER_LIST_SET=()
SIB_AFFECTED_LIST_SET=()

SIBLING_SAMPLE_NAMES=()
SIBLING_GENDERS=()
SIBLING_AFFECTED=()

FAMILY_ID=""
for SAMPLE_ID in ${SAMPLES_LIST[@]}
do
    FAMILY_ID=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].family)"))
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
INPUT_DATA_DIR="${COHORT_WORKFLOW_DIR}/${FAMILY_ID}_pedigree_outstore"
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
    SIB_BAM_LIST+="'${INPUT_DATA_DIR}/${SIBLING_SAMPLE_NAMES[$n]}_merged.bam' "
    SIB_BAI_LIST+="'${INPUT_DATA_DIR}/${SIBLING_SAMPLE_NAMES[$n]}_merged.bam.bai' "
done

# Build trio .ped file for parental graph construction
TRIO_PED_FILE="${PROBAND_SAMPLE_NAME}.trio.ped"
rm -f ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "#Family\tID\tFather\tMother\tSex[1=M]\tAffected[2=A]" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "${FAMILY_ID}\t${PROBAND_SAMPLE_NAME}\t${PATERNAL_SAMPLE_NAME}\t${MATERNAL_SAMPLE_NAME}\t$((${SIB_GENDER_LIST[0]} + 1))\t2" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "${FAMILY_ID}\t${PATERNAL_SAMPLE_NAME}\t0\t0\t1\t1" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "${FAMILY_ID}\t${MATERNAL_SAMPLE_NAME}\t0\t0\t2\t1" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
chmod 2770 "${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}"

if [[ ${COHORT_WORKFLOW_DIR} = *[[:space:]]* ]] || [ -z ${COHORT_WORKFLOW_DIR} ]; then
    echo "ERROR: ${COHORT_WORKFLOW_DIR} argument value contains whitespace or is empty"
    exit 1
fi
if [[ ${FAMILY_ID} = *[[:space:]]* ]] || [ -z ${FAMILY_ID} ]; then
    echo "ERROR: ${FAMILY_ID} argument value contains whitespace or is empty"
    exit 1
fi

if [ ! -d "${COHORT_WORKFLOW_DIR}/${FAMILY_ID}_analysis_outstore" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/${FAMILY_ID}_analysis_outstore"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/${FAMILY_ID}_analysis_outstore"
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
echo "toil-vg generate-config --whole_genome >config.cfg" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
## Adjust the containers in the config file that are dependent on genome reference version
echo "sed -i'' config.cfg -e \"s|^preprocess-mem:.*|preprocess-mem: \'50G\'|g\"" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
if [ $GRCh38_REFERENCE_VERSION == false ]; then
    echo "sed -i'' config.cfg -e \"s|_grch38||g\"" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
else
    echo "sed -i'' config.cfg -e \"s|vcf2shebang-docker: \'quay.io/cmarkello/vcf2shebang:latest\'|vcf2shebang-docker: \'quay.io/cmarkello/vcf2shebang_grch38:latest\'|g\" -e \"s|bmtb-docker: \'quay.io/cmarkello/bmtb:latest\'|bmtb-docker: \'quay.io/cmarkello/bmtb_grch38:latest\'|g\"" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
fi

if [ $RESTART == false ]; then
    echo "toil clean ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_jobstore" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh
fi

RESTART_ARG=""
if [ $RESTART == true ]; then
    RESTART_ARG="--restart"
fi

GENOME_BUILD_ARG=""
if [ $GRCh38_REFERENCE_VERSION == false ]; then
    GENOME_BUILD_ARG="--genome_build 'GRCh37'"
else
    GENOME_BUILD_ARG="--genome_build 'GRCh38'"
fi

echo "toil-vg analysis \\
${RESTART_ARG} \\
${GENOME_BUILD_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 60 \\
--rescueJobsFrequency 60 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir onSuccess \\
--config ${COHORT_WORKFLOW_DIR}/config.cfg \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_jobstore \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_outstore \\
--cohort_vcf ${INPUT_DATA_DIR}/${PROBAND_SAMPLE_NAME}.snpeff.unrolled.vcf.gz \\
--trio_pedfile ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}.trio.ped \\
--sample_name ${PROBAND_SAMPLE_NAME} \\
--maternal_name ${MATERNAL_SAMPLE_NAME} \\
--paternal_name ${PATERNAL_SAMPLE_NAME} \\
--sibling_names ${SIB_ID_LIST[@]} \\
--sibling_genders ${SIB_GENDER_LIST[@]} \\
--sibling_affected ${SIB_AFFECTED_LIST[@]} \\
--maternal_bam ${INPUT_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_merged.bam \\
--maternal_bai ${INPUT_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_merged.bam.bai \\
--paternal_bam ${INPUT_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_merged.bam \\
--paternal_bai ${INPUT_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_merged.bam.bai \\
--siblings_bam ${SIB_BAM_LIST[@]} \\
--siblings_bai ${SIB_BAI_LIST[@]} \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_analysis_workflow.sh

exit
