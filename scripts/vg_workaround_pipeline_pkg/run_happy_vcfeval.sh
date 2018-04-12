#!/bin/bash

## Example run: ./run_happy_vcfeval.sh UDP10618 UDP10618 /data/markellocj/ t146

COHORT_NAME=$1
SAMPLE_NAME=$2
WORK_DIR=$3
VG_VERSION=$4

UDN_VCF_NAME=`ls /data/Udpdata/Individuals/${SAMPLE_NAME}/WGS/Aligned/HudsonAlpha_1/ | grep '^UDN' | cut -f 1 -d '.'`
VCF_ANALYSIS_DIR="${WORK_DIR}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare"
TRUE_SNPCHIP_VCF_FILE="${WORK_DIR}/GSfinal_report/${COHORT_NAME}/${SAMPLE_NAME}_only_SNPchip.vcf.gz"

UDN_BASELINE_FILE="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/${UDN_VCF_NAME}.new.vcf.gz"
SNP_BASELINE_FILE="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/${SAMPLE_NAME}_only_SNPchip.new.vcf.gz"
CALLED_FILE="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/${COHORT_NAME}_dragen_call_from_vg_surject_vg_${VG_VERSION}/${SAMPLE_NAME}.new.vcf.gz"
SDF_DIR="${HOME}/sdf_references/hg19.fa.sdf"
VCFEVAL_UDN_OUTPUT="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/vcfeval_output_UDN"
VCFEVAL_SNP_OUTPUT="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/vcfeval_output_SNP"
VCFEVAL_UDN_vs_SNP_OUTPUT="${HOME}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/output_${SAMPLE_NAME}_vg_${VG_VERSION}/${SAMPLE_NAME}_vcf_compare/vcfeval_output_UDN_vs_SNP"

module load vt/0.5
module load hap.py

mkdir ${VCF_ANALYSIS_DIR}

# Preprocess UDN vcf
cp /data/Udpdata/Individuals/${SAMPLE_NAME}/WGS/Aligned/HudsonAlpha_1/${UDN_VCF_NAME}.vcf.gz ${VCF_ANALYSIS_DIR}
cd ${VCF_ANALYSIS_DIR}
bgzip -d ${UDN_VCF_NAME}.vcf.gz
awk '{if($0 !~ /^#|^chr/) print "chr"$0; else print $0}' ${UDN_VCF_NAME}.vcf > ${UDN_VCF_NAME}.new.vcf
bgzip ${UDN_VCF_NAME}.new.vcf
tabix -f -p vcf ${UDN_VCF_NAME}.new.vcf.gz
bgzip ${UDN_VCF_NAME}.vcf
tabix -f -p vcf ${UDN_VCF_NAME}.vcf.gz

# Preprocess SNP chip vcf
cp /data/markellocj/GSfinal_report/${COHORT_NAME}/${SAMPLE_NAME}_only_SNPchip.vcf.gz ${VCF_ANALYSIS_DIR}
cd ${VCF_ANALYSIS_DIR}
bgzip -d ${SAMPLE_NAME}_only_SNPchip.vcf.gz
awk '{if($0 !~ /^#|^chr/) print "chr"$0; else print $0}' ${SAMPLE_NAME}_only_SNPchip.vcf > ${SAMPLE_NAME}_only_SNPchip.new.vcf
bgzip ${SAMPLE_NAME}_only_SNPchip.new.vcf
tabix -f -p vcf ${SAMPLE_NAME}_only_SNPchip.new.vcf.gz
bgzip ${SAMPLE_NAME}_only_SNPchip.vcf
tabix -f -p vcf ${SAMPLE_NAME}_only_SNPchip.vcf.gz

# Preprocess vg called vcf
cp ${WORK_DIR}/${COHORT_NAME}_run/output_vg_${VG_VERSION}/${COHORT_NAME}_dragen_call_from_vg_surject_vg_${VG_VERSION}/${SAMPLE_NAME}.vcf.gz ${VCF_ANALYSIS_DIR}
cd ${VCF_ANALYSIS_DIR}
bgzip -d ${SAMPLE_NAME}.vcf.gz
awk '{if($0 !~ /^#|^chr/) print "chr"$0; else print $0}' ${SAMPLE_NAME}.vcf > ${SAMPLE_NAME}.new.vcf
bgzip ${SAMPLE_NAME}.new.vcf
tabix -f -p vcf ${SAMPLE_NAME}.new.vcf.gz
bgzip ${SAMPLE_NAME}.vcf
tabix -f -p vcf ${SAMPLE_NAME}.vcf.gz



## Run hap.py against UDN linear vcf
cd ${VCF_ANALYSIS_DIR}
hap.py ${UDN_VCF_NAME}.new.vcf.gz ${SAMPLE_NAME}.new.vcf.gz -o happy_${SAMPLE_NAME}_UDNvcf --threads 32 2> happy_output_UDNvcf.txt

## Run hap.py against SNP chip vcf
cd ${VCF_ANALYSIS_DIR}
hap.py ${SAMPLE_NAME}_only_SNPchip.new.vcf.gz ${SAMPLE_NAME}.new.vcf.gz -o happy_${SAMPLE_NAME}_snpchip --threads 32 2> happy_output_snpchip.txt

## Run hap.py UDN vcf against SNP chip vcf
cd ${VCF_ANALYSIS_DIR}
hap.py ${SAMPLE_NAME}_only_SNPchip.new.vcf.gz ${UDN_VCF_NAME}.new.vcf.gz -o happy_${SAMPLE_NAME}_UDN_snpchip --threads 32 2> happy_output_UDN_snpchip.txt


module load singularity

## Run vcfeval against UDN linear vcf
cd ${WORK_DIR}
rm -fr ${VCF_ANALYSIS_DIR}/vcfeval_output
singularity exec -H ${PWD}:${HOME}  docker://realtimegenomics/rtg-tools rtg vcfeval -T 32 -b ${UDN_BASELINE_FILE} -c ${CALLED_FILE} -t ${SDF_DIR} -o ${VCFEVAL_UDN_OUTPUT}

## Run vcfeval against SNP linear vcf
cd ${WORK_DIR}
rm -fr ${VCFEVAL_SNP_OUTPUT}
singularity exec -H ${PWD}:${HOME}  docker://realtimegenomics/rtg-tools rtg vcfeval -T 32 -b ${SNP_BASELINE_FILE} -c ${CALLED_FILE} -t ${SDF_DIR} -o ${VCFEVAL_SNP_OUTPUT}

## Run vcfeval UDN vcf against SNP chip vcf
cd ${WORK_DIR}
rm -fr ${VCFEVAL_UDN_vs_SNP_OUTPUT}
singularity exec -H ${PWD}:${HOME}  docker://realtimegenomics/rtg-tools rtg vcfeval -T 32 -b ${SNP_BASELINE_FILE} -c ${UDN_BASELINE_FILE} -t ${SDF_DIR} -o ${VCFEVAL_UDN_vs_SNP_OUTPUT}

