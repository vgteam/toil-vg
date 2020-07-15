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
    -m Mother UDP ID (in format UDP####)
    -f Father UDP ID (in format UDP####)
    -s List of Sibling UDP ID, Proband ID must be first in the list (in format UDP#### UDP#### UDP#### ...)
    -c PATH to .ped file containing only the mother-father-proband trio samples
    -w PATH to where the UDP cohort will be processed and where the input reads will be stored
    -g PATH to the workflow input directory
    -i List of Sibling Gender IDs. 0=male, 1=female. Must be same order as the input to -s argument.
    -b List of Sibling affected status. 0=unaffected, 1=affected. Must be same order as the input to -s argument.
    -a PATH to chromosome annotation directory used by vcftoshebang.
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
while getopts "m:f:s:c:w:g:i:b:a:e:d:v:r:t:h" OPTION; do
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
        c)
            TRIO_PEDIGREE_FILE=$OPTARG
        ;;
        w)
            COHORT_WORKFLOW_DIR=$OPTARG
        ;;
        g)
            WORKFLOW_INPUT_DIR=$OPTARG
        ;;
        i)
            SIBLING_GENDERS+=($OPTARG)
        ;;
        b)
            SIBLING_AFFECTED+=($OPTARG)
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
READ_DATA_DIR="${COHORT_WORKFLOW_DIR}/input_reads"
SIB_READ_PAIR_LIST=""
SIB_ID_LIST=""
SIB_GENDER_LIST=""
SIB_AFFECT_LIST=""
SIBLING_SAMPLE_NAMES_LEN=${#SIBLING_SAMPLE_NAMES[@]}
if [ ${#SIBLING_SAMPLE_NAMES[@]} -gt 1 ]; then
    for SIBLING_ID in ${SIBLING_SAMPLE_NAMES[@]:1}
    do
      SIB_READ_PAIR_LIST+="${READ_DATA_DIR}/${SIBLING_ID}_read_pair_1.fq.gz ${READ_DATA_DIR}/${SIBLING_ID}_read_pair_2.fq.gz "
      SIB_ID_LIST+="'${SIBLING_ID}' "
    done
fi

for (( n=1; n<=${#SIBLING_SAMPLE_NAMES[@]}; n++ ))
do
    SIB_GENDER_LIST+="'${SIBLING_GENDERS[$n]}' "
    SIB_AFFECT_LIST+="'${SIBLING_AFFECTED[$n]}' "
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
    rm -fr ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore
    rm -fr ${COHORT_WORKFLOW_DIR}/tmp
fi

if [ ! -d "${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore"
fi
if [ ! -d "${COHORT_WORKFLOW_DIR}/tmp" ]; then
    mkdir -p "${COHORT_WORKFLOW_DIR}/tmp"
    chmod 2770 "${COHORT_WORKFLOW_DIR}/tmp"
fi

if [ ! -d "/data/$USER/singularity_cache" ]; then
    mkdir -p "/data/$USER/singularity_cache"
    chmod 2770 "/data/$USER/singularity_cache"
fi

rm -f ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo '#!/bin/bash' >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "module load singularity python/3.7" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "source ${TOIL_VG_DIR}/toilvg_venv/bin/activate" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "export TOIL_SLURM_ARGS='-t 20:00:00'" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "export SINGULARITY_CACHEDIR=/data/$USER/singularity_cache" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
echo "cd ${COHORT_WORKFLOW_DIR}" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
if [ $RESTART == false ]; then
    echo "toil clean ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_jobstore" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
fi
RESTART_ARG=""
if [ $RESTART == true ]; then
    RESTART_ARG="--restart"
fi
if [ $RUN_SMALL_TEST == false ]; then
    echo "toil-vg pedigree \\
${RESTART_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 30 \\
--rescueJobsFrequency 30 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir onSuccess \\
--whole_genome_config \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_jobstore \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore \\
${PROBAND_SAMPLE_NAME} \\
${MATERNAL_SAMPLE_NAME} \\
${PATERNAL_SAMPLE_NAME} \\
--sibling_names ${SIB_ID_LIST[@]} \\
--sibling_genders ${SIB_GENDER_LIST[@]} \\
--sibling_affected ${SIB_AFFECT_LIST[@]} \\
--fastq_proband ${READ_DATA_DIR}/${PROBAND_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PROBAND_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_maternal ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_paternal ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_siblings ${SIB_READ_PAIR_LIST[@]} \\
--reads_per_chunk 20000000 \\
--ref_fasta ${WORKFLOW_INPUT_DIR}/hs37d5.fa \\
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai \\
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/hs37d5.dict \\
--xg_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.xg \\
--gcsa_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gcsa \\
--gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_decoys.gbwt \\
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \\
--path_list ${WORKFLOW_INPUT_DIR}/path_list_whole_genome.txt \\
--ped_file ${TRIO_PEDIGREE_FILE} \\
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip \\
--genetic_map ${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar \\
--bam_output \\
--use_decoys \\
--force_phasing True \\
--indel_realign_bams \\
--snpeff_annotation \\
--run_dragen \\
--dragen_ref_index_name 'hs37d5_v7' \\
--udp_data_dir 'Udpbinfo' \\
--run_analysis \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR} \\
--helix_username $USER" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
else
    echo "toil-vg pedigree \\
${RESTART_ARG} \\
--setEnv PATH=\$PATH \\
--batchSystem Slurm \\
--statePollingWait 30 \\
--rescueJobsFrequency 30 \\
--container Singularity \\
--logInfo \\
--logFile ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.log \\
--workDir ${COHORT_WORKFLOW_DIR}/tmp \\
--cleanWorkDir always \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_jobstore \\
${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_outstore \\
${PROBAND_SAMPLE_NAME} \\
${MATERNAL_SAMPLE_NAME} \\
${PATERNAL_SAMPLE_NAME} \\
--sibling_names ${SIB_ID_LIST[@]} \\
--sibling_genders ${SIB_GENDER_LIST[@]} \\
--sibling_affected ${SIB_AFFECT_LIST[@]} \\
--fastq_proband ${READ_DATA_DIR}/${PROBAND_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PROBAND_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_maternal ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${MATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_paternal ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_1.fq.gz ${READ_DATA_DIR}/${PATERNAL_SAMPLE_NAME}_read_pair_2.fq.gz \\
--fastq_siblings ${SIB_READ_PAIR_LIST[@]} \\
--ref_fasta ${WORKFLOW_INPUT_DIR}/hs37d5.fa \\
--ref_fasta_index ${WORKFLOW_INPUT_DIR}/hs37d5.fa.fai \\
--ref_fasta_dict ${WORKFLOW_INPUT_DIR}/hs37d5.dict \\
--xg_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.xg \\
--gcsa_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gcsa \\
--gbwt_index ${WORKFLOW_INPUT_DIR}/snp1kg_maf0.01_chr21.gbwt \\
--id_ranges ${WORKFLOW_INPUT_DIR}/path_list_21.txt \\
--path_list ${WORKFLOW_INPUT_DIR}/path_list_21.txt \\
--ped_file ${TRIO_PEDIGREE_FILE} \\
--snpeff_database ${WORKFLOW_INPUT_DIR}/snpEff_v4_3_GRCh37.75.zip \\
--genetic_map ${WORKFLOW_INPUT_DIR}/genetic_map_GRCh37.tar \\
--bam_output \\
--force_phasing True \\
--indel_realign_bams \\
--snpeff_annotation \\
--run_dragen \\
--dragen_ref_index_name 'hs37d5_v7' \\
--udp_data_dir 'Udpbinfo' \\
--run_analysis \\
--chrom_dir ${CHROM_ANNOT_DIR} \\
--edit_dir ${EDIT_ANNOT_DIR} \\
--cadd_data ${CADD_DATA_DIR} \\
--helix_username $USER" >> ${COHORT_WORKFLOW_DIR}/${PROBAND_SAMPLE_NAME}_pedigree_workflow.sh
fi

exit

