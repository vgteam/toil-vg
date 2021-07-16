#!/bin/bash
#################################################################################################
##
##  Script to setup a bash script to run the entire toil-vg pedigree pipeline.
##
##  Inputs:
##
##  Assumptions:
##      The UDP cohort ID is the same as the UDP ID for the proband sample.
##
##  Last modified:
##  Last modified by: Charles Markello
##
#################################################################################################

## Create help statement
usage(){
cat << EOF

This script setups up a bash script to run a UDP cohort through all stages of the toil-vg
pedigree pipeline on the NIH Biowulf Cluster.

Inputs:
Inputs:
    -f Cohort name. Should be the same as the proband sample name.
    -c PATH to .ped file containing all samples in the family.
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -a PATH to chromosome annotation directory used by vcftoshebang.
    -e PATH to directory containing master edit files used by vcftoshebang.
    -d PATH to cadd engine data directory.
    -v PATH to the toil_vg repository
    -i (OPTIONAL, default=false) Set to 'true' to run workflow using the udp nih illumina dragen module for variant calling
    -r (OPTIONAL, default=false) Set to 'true' to restart an incompletely ran workflow
    -b (OPTIONAL, default=false) Set to 'true' to setup a workflow that runs GRCh38-based input
    -t (OPTIONAL, default=false) Set to 'true' if running workflow on small HG002 chr21 test data
    
Outputs:

Assumptions:

EOF

}

## Check number of arguments
if [ $# -lt 8 ] || [[ $@ != -* ]]; then
    usage
    exit 1
fi

## DEFAULT PARAMETERS
USE_DRAGEN=false
RUN_SMALL_TEST=false
RESTART=false
GRCh38_REFERENCE_VERSION=false

## Parse through arguments
while getopts "f:c:w:g:a:e:d:v:i:r:b:t:h" OPTION; do
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
        g)
            WORKFLOW_INPUT_DIR=$OPTARG
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
        i)
            USE_DRAGEN=$OPTARG
        ;;
        r)
            RESTART=$OPTARG
        ;;
        b)
            GRCh38_REFERENCE_VERSION=$OPTARG
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

# Extract sample information from input family .ped file
source ${TOIL_VG_DIR}/toilvg_venv/bin/activate
pip install ped_parser
TRIO_PED_FILE="${COHORT_NAME}.trio.ped"

READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
SAMPLES_LIST=($(python3 -c "import ped_parser; print(list(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals.keys()))" | tr -d '[],' | tr -d \'))
PROBAND_NAME=""
OFFSPRING_GENDER_LIST=()
OFFSPRING_AFFECTED_LIST=()
SIB_READ_PAIR_LIST=""
SIB_ID_LIST=""
SIB_ID_LIST_SET=()
SIB_GENDER_LIST_SET=()
SIB_AFFECTED_LIST_SET=()

for SAMPLE_ID in ${SAMPLES_LIST[@]}
do
    SAMPLE_MOM=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].mother)"))
    SAMPLE_DAD=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].father)"))
    SAMPLE_GENDER=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].sex)"))
    SAMPLE_AFFECTED=($(python3 -c "import ped_parser; print(ped_parser.FamilyParser(family_info=open('${COHORT_PED_FILE}','r'), family_type='ped').individuals['${SAMPLE_ID}'].phenotype)"))
    if [[ "$SAMPLE_ID" == "$COHORT_NAME" ]]; then
        PROBAND_NAME=${SAMPLE_ID}
        OFFSPRING_GENDER_LIST+=($((${SAMPLE_GENDER} - 1)))
        if [[ ${SAMPLE_AFFECTED} -eq 1 ]]; then
            OFFSPRING_AFFECTED_LIST+=(1)
        else
            OFFSPRING_AFFECTED_LIST+=(0)
        fi
        MATERNAL_SAMPLE_NAME="${SAMPLE_MOM}"
        PATERNAL_SAMPLE_NAME="${SAMPLE_DAD}"
    elif [[ ${SAMPLE_MOM} != '0' ]]; then
        SIB_ID_LIST_SET+=(${SAMPLE_ID})
        SIB_READ_PAIR_LIST+="${READ_DATA_DIR}/${SAMPLE_ID}_read_pair_1.fq.gz ${READ_DATA_DIR}/${SAMPLE_ID}_read_pair_2.fq.gz "
        SIB_GENDER_LIST+=(${SAMPLE_GENDER})
        SIB_AFFECTED_LIST+=(${SAMPLE_AFFECTED})
    fi
done

for (( n=0; n<${#SIB_ID_LIST_SET[@]}; n++ ))
do
    OFFSPRING_GENDER_LIST+=($((${SIB_GENDER_LIST[$n]} - 1)))
    if [[ ${SIB_AFFECTED_LIST[$n]} -eq 1 ]]; then
        OFFSPRING_AFFECTED_LIST+=(1)
    else
        OFFSPRING_AFFECTED_LIST+=(0)
    fi
done

if [ ${#SIB_ID_LIST_SET[@]} -gt 0 ]; then
    SIB_READ_PAIR_LIST="--fastq_siblings "
    SIB_ID_LIST="--sibling_names "
    for SIBLING_ID in ${SIB_ID_LIST_SET[@]}
    do
      SIB_READ_PAIR_LIST+="${READ_DATA_DIR}/${SIBLING_ID}_read_pair_1.fq.gz ${READ_DATA_DIR}/${SIBLING_ID}_read_pair_2.fq.gz "
      SIB_ID_LIST+="${SIBLING_ID} "
    done
fi

# Make sure directory paths contain no white-space
if [[ ${COHORT_WORKFLOW_DIR} = *[[:space:]]* ]]; then
    echo "ERROR: ${COHORT_WORKFLOW_DIR} argument value contains whitespace"
    exit 1
fi
if [[ ${COHORT_NAME} = *[[:space:]]* ]]; then
    echo "ERROR: ${COHORT_NAME} argument value contains whitespace"
    exit 1
fi

# Build trio .ped file for parental graph construction
rm -f ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "#Family\tID\tFather\tMother\tSex[1=M]\tAffected[2=A]" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "${PROBAND_NAME}\t${PROBAND_NAME}\t${PATERNAL_SAMPLE_NAME}\t${MATERNAL_SAMPLE_NAME}\t$((${OFFSPRING_GENDER_LIST[0]} + 1))\t2" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "${PROBAND_NAME}\t${PATERNAL_SAMPLE_NAME}\t0\t0\t1\t1" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
echo -e "${PROBAND_NAME}\t${MATERNAL_SAMPLE_NAME}\t0\t0\t2\t1" >> ${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}
chmod 2770 "${COHORT_WORKFLOW_DIR}/${TRIO_PED_FILE}"
deactivate

# Make outstore, jobstore, and container cache directories
if [ ! -d "${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_outstore" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_outstore"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_outstore"
fi
if [ ! -d "${COHORT_WORKFLOW_DIR}/tmp" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/tmp"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/tmp"
fi

if [ ! -d "/data/$USER/singularity_cache" ]; then
    mkdir -p "/data/$USER/singularity_cache"
    chmod 2770 "/data/$USER/singularity_cache"
fi

# Build the workflow script
rm -f ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
echo "module load singularity python/3.7" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
echo "source ${TOIL_VG_DIR}/toilvg_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
echo "export TOIL_SLURM_ARGS='-t 20:00:00'" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
echo "export SINGULARITY_CACHEDIR=/data/$USER/singularity_cache" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
if [ $RESTART == false ]; then
    echo "toil clean ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_jobstore" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
fi
RESTART_ARG=""
if [ $RESTART == true ]; then
    RESTART_ARG="--restart"
fi

CALLER="deepvariant"
DRAGEN_ARGS=""
if [ $USE_DRAGEN == true ]; then
    CALLER="dragen"
    DRAGEN_ARGS="--dragen_ref_index_name 'hs37d5_v7' --udp_data_dir 'Udpbinfo'"
fi

if [ $GRCh38_REFERENCE_VERSION == false ]; then
    if [ $RUN_SMALL_TEST == false ]; then
    echo "toil-vg pedigree \\
${RESTART_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 120 \\
--rescueJobsFrequency 120 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir onSuccess \\
--whole_genome_config \\
${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_jobstore \\
${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_outstore \\
${PROBAND_NAME} \\
${MATERNAL_SAMPLE_NAME} \\
${PATERNAL_SAMPLE_NAME} \\
${SIB_ID_LIST} \\
--sibling_genders ${OFFSPRING_GENDER_LIST[@]} \\
--sibling_affected ${OFFSPRING_AFFECTED_LIST[@]} \\
--fastq_proband ${READ_DATA_DIR}/${PROBAND_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PROBAND_NAME}_read_pair_2.fq.gz \\
--fastq_maternal ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_paternal ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
${SIB_READ_PAIR_LIST} \\
--reads_per_chunk 200000000 \\
--ref_fasta ${WORKFLOW_INPUT_DIR}/hs37d5.fa \\
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai \\
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/hs37d5.dict \\
--caller ${CALLER} \\
--mapper giraffe \\
--use_haplotypes \\
--xg_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.xg \\
--gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gbwt \\
--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gg \\
--minimizer_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.min \\
--distance_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.dist \\
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \\
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \\
--ped_file ${TRIO_PED_FILE} \\
--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data.tar.gz \\
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh37.75 \\
--genetic_map ${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar \\
--bam_output \\
--use_decoys \\
--indel_realign_bams \\
--snpeff_annotation \\
${DRAGEN_ARGS} \\
--run_analysis \\
--cadd_lines 100000 \\
--split_lines 100000 \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR} \\
--helix_username $USER" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
    else
    echo "toil-vg pedigree \\
${RESTART_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 120 \\
--rescueJobsFrequency 120 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir always \\
${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_jobstore \\
${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_outstore \\
${PROBAND_NAME} \\
${MATERNAL_SAMPLE_NAME} \\
${PATERNAL_SAMPLE_NAME} \\
${SIB_ID_LIST} \\
--sibling_genders ${OFFSPRING_GENDER_LIST[@]} \\
--sibling_affected ${OFFSPRING_AFFECTED_LIST[@]} \\
--fastq_proband ${READ_DATA_DIR}/${PROBAND_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PROBAND_NAME}_read_pair_2.fq.gz \\
--fastq_maternal ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_paternal ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
${SIB_READ_PAIR_LIST} \\
--ref_fasta ${WORKFLOW_INPUT_DIR}/hs37d5.fa \\
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai \\
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/hs37d5.dict \\
--xg_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.xg \\
--gcsa_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gcsa \\
--gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gbwt \\
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_21.txt \\
--path_list ${WORKFLOW_INPUT_DIR}/path_list_21.txt \\
--ped_file ${TRIO_PED_FILE} \\
--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data.tar.gz \\
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh37.75 \\
--genetic_map ${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar \\
--bam_output \\
--force_phasing True \\
--indel_realign_bams \\
--snpeff_annotation \\
${DRAGEN_ARGS} \\
--run_analysis \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR} \\
--helix_username $USER" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
    fi
else
    echo "toil-vg pedigree \\
${RESTART_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 120 \\
--rescueJobsFrequency 120 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir onSuccess \\
--whole_genome_config \\
${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_jobstore \\
${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_outstore \\
${PROBAND_NAME} \\
${MATERNAL_SAMPLE_NAME} \\
${PATERNAL_SAMPLE_NAME} \\
${SIB_ID_LIST} \\
--sibling_genders ${OFFSPRING_GENDER_LIST[@]} \\
--sibling_affected ${OFFSPRING_AFFECTED_LIST[@]} \\
--fastq_proband ${READ_DATA_DIR}/${PROBAND_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PROBAND_NAME}_read_pair_2.fq.gz \\
--fastq_maternal ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_paternal ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
${SIB_READ_PAIR_LIST} \\
--reads_per_chunk 200000000 \\
--ref_fasta ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna \\
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.fna.fai \\
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.compact_decoys.dict \\
--caller ${CALLER} \\
--mapper giraffe \\
--use_haplotypes \\
--xg_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.xg \\
--gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gbwt \\
--graph_gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.gg \\
--minimizer_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.min \\
--distance_index ${WORKFLOW_INPUT_DIR}/snp1kg_decoys.dist \\
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \\
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \\
--ped_file ${TRIO_PED_FILE} \\
--eagle_data ${WORKFLOW_INPUT_DIR}/eagle_data_grch38.tar.gz \\
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v5_0_GRCh38.99.zip \\
--bam_output \\
--use_decoys \\
--indel_realign_bams \\
--snpeff_annotation \\
${DRAGEN_ARGS} \\
--run_analysis \\
--cadd_lines 100000 \\
--split_lines 100000 \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR} \\
--helix_username $USER" >> ${COHORT_WORKFLOW_DIR}/${COHORT_NAME}_pedigree_workflow.sh
fi

exit

