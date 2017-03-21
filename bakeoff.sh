#!/usr/bin/env bash

## TODO: Get this going in some sort of CI framework
##       Ideally:  vg merge to master triggers docker build then a run of this script
##                 with the accuracy results logged to some kind of persistent leaderboard
##
## For now, all we have is the basic logic to do one run of bakeoff, expecting
## virtualenv etc has already been set up. 

usage() { printf "Usage: $0 [Options] <Output-prefix> <Ouptut F1 File.tsv> \nOptions:\n\t-m\tmesos\n\t-f\tfast (just BRCA1)\n\t-D\tno Docker\n\t-t <N>\tthreads [DEFAULT=7]\n\t-c <F>\tconfig file\n\t-e <EVAL-TSV>\t run mapeval\n\n" 1>&2; exit 1; }

FAST=0
MESOS=0
DOCKER=1
CORES=7
CONFIG=0
while getopts "fmDt:c:e:" o; do
    case "${o}" in
        f)
            FAST=1
            ;;
        m)
            MESOS=1
            ;;
        D)
            DOCKER=0
            ;;
        t)
            CORES=$OPTARG
            ;;
        c)
            CONFIG=$OPTARG
            ;;
        e)
            MEFILE=${OPTARG}
            ;;		  
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [ "$#" -ne 2 ]; then
	 usage
	 exit 1
fi

# all input data expected to be here:
BAKEOFF_STORE="s3://glennhickey-bakeoff-store"

PREFIX=$1
F1FILE=$2

# General Options
GEN_OPTS="--realTimeLogging --logInfo"
if [ "$CONFIG" != "0" ]; then
	 GEN_OPTS="${GEN_OPTS} --config ${CONFIG}"
fi
INDEX_OPTS="--gcsa_index_cores ${CORES} --kmers_cores ${CORES}"
MAP_OPTS="--alignment_cores ${CORES} --interleaved"
MAP_OPTS_SE="--alignment_cores ${CORES}"
CALL_OPTS="--calling_cores ${CORES}"
VCFEVAL_OPTS="--vcfeval_cores ${CORES}"
RUN_OPTS="${INDEX_OPTS} ${MAP_OPTS} ${CALL_OPTS} ${VCFEVAL_OPTS}"
SIM_OPTS="--sim_chunks ${CORES} --seed 0 --gam"
NUM_SIM_READS=100000

# Hack in support for switching between mesos and local here
# (Note, for job store and out store, we will tack on -REGION to make them unique)
if [ "$MESOS" == "1" ]; then
	 JOB_STORE="aws:us-west-2:${PREFIX}-bakeoff-job-store"
	 OUT_STORE="aws:us-west-2:${PREFIX}-bakeoff-out-store"
	 CP_OUT_STORE="s3://${PREFIX}-bakeoff-out-store"
	 BS_OPTS="--batchSystem=mesos --mesosMaster=mesos-master:5050"
	 CP_CMD="aws s3 cp"
	 RM_CMD="aws s3 rm"
else
	 JOB_STORE="./${PREFIX}-bakeoff-job-store"
	 OUT_STORE="./${PREFIX}-bakeoff-out-store"
	 CP_OUT_STORE=${OUT_STORE}
	 BS_OPTS=""
	 CP_CMD="cp"
	 RM_CMD="rm"
fi

# Hack in support for turning Docker on and off
if [ "$DOCKER" == "0" ]; then
    DOCKER_OPTS="--container None"
else
    DOCKER_OPTS=""
fi

# Run toil-vg on a bakeoff region
function run_bakeoff_region {
	 local REGION=$1
	 local CHROM=$2
	 local OFFSET=$3

	 # erase the job store and f1 output
	 toil clean ${JOB_STORE}-${REGION,,}
	 F1_SCORE="${CP_OUT_STORE}-${REGION,,}/NA12878_vcfeval_output_f1.txt"
	 $RM_CMD $F1_SCORE

	 # run the whole pipeline
	 toil-vg run ${JOB_STORE}-${REGION,,} NA12878 ${OUT_STORE}-${REGION,,} --fastq ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.fq.gz --chroms ${CHROM} --call_opts "--offset ${OFFSET}" --graphs ${BAKEOFF_STORE}/snp1kg-${REGION}.vg --vcfeval_baseline ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.vcf.gz --vcfeval_fasta ${BAKEOFF_STORE}/chrom.fa.gz ${BS_OPTS} ${DOCKER_OPTS} ${GEN_OPTS} ${RUN_OPTS}

	 toil clean ${JOB_STORE}-${REGION,,}

	 # harvest the f1 output and append it to our table
	 rm -f temp_f1.txt
	 $CP_CMD $F1_SCORE ./temp_f1.txt
	 if [ -f "./temp_f1.txt" ] ; then
		  printf "${REGION}\t$(cat ./temp_f1.txt)\n" >> ${F1FILE}
	 else
		  printf "${REGION}\tERROR\n" >> ${F1FILE}
	 fi

	 # optionally run the simulation benchmark
	 if test ${MEFILE+defined}; then
		  run_mapeval $REGION $CHROM
	 fi
}

# Run simulation mapping evaluation
function run_mapeval {
	 local REGION=$1
	 local CHROM=$2

	 # erase the job store and f1 output
	 toil clean ${JOB_STORE}-${REGION,,}
	 STATS_TSV="${CP_OUT_STORE}-me-${REGION,,}/stats.tsv"
	 $RM_CMD $STATS_TSV

	 # simulate some reads
	 toil-vg sim ${JOB_STORE}-${REGION,,} ${OUT_STORE}-${REGION,,}/genome.xg ${NUM_SIM_READS} ${OUT_STORE}-me-${REGION,,} ${BS_OPTS} ${DOCKER_OPTS} ${SIM_OPTS}
	 toil clean ${JOB_STORE}-${REGION,,}

	 # generate mapeval results
	 toil-vg mapeval ${JOB_STORE}-${REGION,,} ${OUT_STORE}-me-${REGION,,}  ${OUT_STORE}-me-${REGION,,}/true.pos ${BS_OPTS} ${BS_OPTS} ${DOCKER_OPTS} ${GEN_OPTS} ${MAP_OPTS_SE} --bwa --bwa-paired --index-bases ${OUT_STORE}-${REGION,,}/genome --gam-names vg --gam_input_reads ${OUT_STORE}-me-${REGION,,}/sim.gam --fasta ${BAKEOFF_STORE}/${REGION}.fa --vg-paired

	 toil clean ${JOB_STORE}-${REGION,,}
	 
	 # harvest the stats output and append it to our table
	 rm -f temp_stats.tsv
	 $CP_CMD $STATS_TSV ./temp_stats.tsv
	 if [ -f "./temp_stats.tsv" ] ; then
		  printf "${REGION}\n$(cat ./temp_stats.tsv)\n" >> ${MEFILE}
	 else
		  printf "${REGION}\nERROR\n" >> ${MEFILE}
	 fi	 
}	 

# Run the regions in series

rm -f $F1FILE
if test ${MEFILE+defined}; then
	 rm -rf $MEFILE
fi

run_bakeoff_region BRCA1 17 43044293
if [ "$FAST" != "1" ]; then
	 run_bakeoff_region BRCA2 13 32314860
	 run_bakeoff_region SMA 5 69216818
	 run_bakeoff_region LRC-KIR 19 54025633
	 run_bakeoff_region MHC 6 28510119
fi

echo "Created the following:"
echo ${JOB_STORE}-brca1
echo ${OUT_STORE}-brca1
if test ${MEFILE+defined}; then
	 echo ${OUT_STORE}-brca1-me
fi
if [ "$FAST" != "1" ]; then
	 echo ${JOB_STORE}-brca2
	 echo ${OUT_STORE}-brca2
	 echo ${JOB_STORE}-sma
	 echo ${OUT_STORE}-sma
	 echo ${JOB_STORE}-lrc-kir
	 echo ${OUT_STORE}-lrc-kir
	 echo ${JOB_STORE}-mhc
	 echo ${OUT_STORE}-mhc
	 if test ${MEFILE+defined}; then
		  echo ${OUT_STORE}-brca2-me
		  echo ${OUT_STORE}-sma-me
		  echo ${OUT_STORE}-lrc-kir-me
		  echo ${OUT_STORE}-mhc-me
	 fi
fi

        
