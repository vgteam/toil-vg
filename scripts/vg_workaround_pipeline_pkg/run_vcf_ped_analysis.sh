#!/bin/bash

COHORT_NAME="UDN622150"
COHORT_WORKDIR="${HOME}/${COHORT_NAME}_run"
JOINT_VCF="${COHORT_WORKDIR}/${COHORT_NAME}_dragen_call_from_vg_surject_vg_t146/joint_genotyped_${COHORT_NAME}.vcf.gz"
OUTPUT="${COHORT_WORKDIR}/mendelian_analysis_vg_t146/${COHORT_NAME}_mendelian_inconsistent.vcf.gz"
PED_FILE="${HOME}/${COHORT_NAME}_run/mendelian_analysis_vg_t146/${COHORT_NAME}.ped"
WORK_DIR="/data/markellocj/"
SDF_DIR="${HOME}/sdf_references/hs37d5.fa.sdf"

module load vt/0.5
module load hap.py
module load singularity

## Run pedigree analysis on vg surjected and Dragen called variants
cd ${WORK_DIR}
singularity exec -H ${PWD}:${HOME}  docker://realtimegenomics/rtg-tools rtg mendelian -i ${JOINT_VCF} --output-inconsistent ${OUTPUT} --pedigree ${PED_FILE} -l -t ${SDF_DIR}


