#!/usr/bin/env bash

## TODO: Get this going in some sort of CI framework
##       Ideally:  vg merge to master triggers docker build then a run of this script
##                 with the accuracy results logged to some kind of persistent leaderboard
##
## For now, all we have is the basic logic to do one run of bakeoff, expecting
## virtualenv etc has already been set up. 

usage() { printf "Usage: $0 [Options] <Output-prefix> \nOptions:\n\t-m\tmesos\n\t-f\tfast (just BRCA1)\n\\t-D\tno Docker\n\t-v <I>\tvg docker image\n\t-t <N>\tthreads [DEFAULT=7]\n\t-c <F>\tconfig file\n\t-e\trun mapeval\n\t-s\twrite toil stats\n" 1>&2; exit 1; }

FAST=0
MESOS=0
DOCKER=1
CORES=8
CONFIG=0
ME=0
STATS=0
VG_DOCKER=0
while getopts "fjmDv:t:c:es" o; do
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
        v)
            VG_DOCKER=$OPTARG
            ;;
        t)
            CORES=$OPTARG
            ;;
        c)
            CONFIG=$OPTARG
            ;;
        e)
            ME=1
            ;;
        s)
            STATS=1
            ;;
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [ "$#" -ne 1 ]; then
	 usage
	 exit 1
fi

# all input data expected to be here:
BAKEOFF_STORE="s3://glennhickey2-bakeoff-store"

# output summary will be here
PREFIX=$1
F1FILE="${PREFIX}-f1.tsv"
if [ "$ME" != 0 ]; then
	 MEFILE="${PREFIX}-mapeval.tsv"
fi
if [ "$STATS" != 0 ]; then
	 STATSFILE="${PREFIX}-stats.tsv"
fi
MAPTIMESFILE="${PREFIX}-maptimes.tsv"
rm -f $MAPTIMESFILE

# General Options
GEN_OPTS="--realTimeLogging --logInfo"
if [ "$CONFIG" != "0" ]; then
	 GEN_OPTS="${GEN_OPTS} --config ${CONFIG}"
fi
if test ${STATSFILE+defined}; then
	 GEN_OPTS="${GEN_OPTS} --stats"
	 rm -f ${STATSFILE}
fi

INDEX_OPTS="--gcsa_index_cores ${CORES} --kmers_cores ${CORES} --prune_opts '-p -l 16 -S -e 3'"
MAP_OPTS="--alignment_cores ${CORES} --interleaved"
MAP_OPTS_SE="--alignment_cores ${CORES}"
CALL_OPTS="--calling_cores ${CORES} --call_opts \"--min-mad 0\""
VCFEVAL_OPTS="--vcfeval_cores ${CORES} --vcfeval_opts \" --ref-overlap\""
RUN_OPTS="${INDEX_OPTS} ${MAP_OPTS} ${CALL_OPTS} ${VCFEVAL_OPTS}"
SIM_OPTS="--sim_chunks ${CORES} --seed 0 --gam"
PRUNE_OPTS=
NUM_SIM_READS=50000

# Hack in support for switching between mesos and local here
# (Note, for job store and out store, we will tack on -REGION to make them unique)
if [ "$MESOS" == "1" ]; then
	 JOB_STORE="aws:us-west-2:${PREFIX}-bakeoff-job-store"
	 OUT_STORE="aws:us-west-2:${PREFIX}-bakeoff-out-store"
	 CP_OUT_STORE="s3://${PREFIX}-bakeoff-out-store"
	 PRIVATE_IP=`ifconfig eth0 |grep "inet addr" |awk '{print $2}' |awk -F: '{print $2}'`
	 # To use spot notes, add: --defaultPreemptable --preemptableNodeType c3.8xlarge:0.85
	 BS_OPTS="--batchSystem=mesos --mesosMaster=${PRIVATE_IP}:5050 --nodeType c3.8xlarge --retry 4 --provisioner aws --minNodes 0 --maxNodes 2"
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
elif [ "$VG_DOCKER" != "0" ]; then
	 DOCKER_OPTS="--vg_docker ${VG_DOCKER}"
else
    DOCKER_OPTS=""
fi

# Run toil-vg on a bakeoff region
function run_bakeoff_region {
	 local REGION=$1
	 local CHROM=$2
	 local OFFSET=$3
	 local GRAPH=$4

	 # erase the job store and f1 output
	 toil clean ${JOB_STORE}-${REGION,,}
	 F1_SCORE="${CP_OUT_STORE}-${REGION,,}-${GRAPH}/NA12878_vcfeval_output_f1.txt"
	 $RM_CMD $F1_SCORE
	 LOGFILE="${PREFIX}-${REGION,,}-${GRAPH}.log"

	 # run the whole pipeline
	 eval toil-vg run ${JOB_STORE}-${REGION,,} NA12878 ${OUT_STORE}-${REGION,,}-${GRAPH} --fastq ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.fq.gz --chroms ${CHROM} --vcf_offsets ${OFFSET} --graphs ${BAKEOFF_STORE}/${GRAPH}-${REGION}.vg --vcfeval_baseline ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.vcf.gz --vcfeval_fasta ${BAKEOFF_STORE}/chr${CHROM}.fa.gz ${BS_OPTS} ${DOCKER_OPTS} ${GEN_OPTS} ${RUN_OPTS} --logFile ${LOGFILE}
	 
	 # extract vg map total time
	 MAPTIME=`extract_map_time ${LOGFILE}`
	 printf "${REGION}\t${GRAPH}\t${MAPTIME}\n" >> $MAPTIMESFILE

	 # harvest stats
	 if test ${STATSFILE+defined}; then
		  echo "STATS FOR ${JOB_STORE}-${REGION,,}" >> ${STATSFILE}
		  toil stats ${JOB_STORE}-${REGION,,} >> ${STATSFILE}
		  printf "\n\n" >> ${STATSFILE}
	 fi

	 #toil clean ${JOB_STORE}-${REGION,,}

	 # harvest the f1 output and append it to our table
	 rm -f temp_f1.txt
	 $CP_CMD $F1_SCORE ./temp_f1.txt
	 if [ -f "./temp_f1.txt" ] ; then
		  printf "${REGION}\t${GRAPH}\t$(cat ./temp_f1.txt)\n" >> ${F1FILE}
	 else
		  printf "${REGION}\t${GRAPH}\tERROR\n" >> ${F1FILE}
	 fi

	 # optionally run the simulation benchmark
	 if test ${MEFILE+defined}; then
		  run_mapeval $REGION $GRAPH ${OUT_STORE}-${REGION,,}-${GRAPH} genome
	 fi
}

# Run simulation mapping evaluation
function run_mapeval {
	 local REGION=$1
	 local GRAPH=$2
	 local IN_STORE=$3
	 local IN_IDX=$4

	 # erase the job store and f1 output
	 toil clean ${JOB_STORE}-${REGION,,}
	 STATS_TSV="${CP_OUT_STORE}-me-${REGION,,}-${GRAPH}/stats.tsv"
	 $RM_CMD $STATS_TSV
	 LOGFILE="${PREFIX}-${REGION,,}-${GRAPH}mapeval.log"

	 # simulate some reads
	 toil-vg sim ${JOB_STORE}-${REGION,,} ${IN_STORE}/${IN_IDX}.xg ${NUM_SIM_READS} ${OUT_STORE}-me-${REGION,,}-${GRAPH} ${BS_OPTS} ${DOCKER_OPTS} ${SIM_OPTS}
	 toil clean ${JOB_STORE}-${REGION,,}

	 # generate mapeval results
	 toil-vg mapeval ${JOB_STORE}-${REGION,,} ${OUT_STORE}-me-${REGION,,}-${GRAPH}  ${OUT_STORE}-me-${REGION,,}-${GRAPH}/true.pos ${BS_OPTS} ${BS_OPTS} ${DOCKER_OPTS} ${GEN_OPTS} ${MAP_OPTS_SE} --bwa --bwa-paired --index-bases ${IN_STORE}/${IN_IDX} --gam-names vg --gam_input_reads ${OUT_STORE}-me-${REGION,,}-${GRAPH}/sim.gam --fasta ${BAKEOFF_STORE}/${REGION}.fa --vg-paired --logFile ${LOGFILE}

	 # extract vg map total time
	 MAPTIME=`extract_map_time ${LOGFILE}`
	 printf "SIM-${REGION}\t${GRAPH}\t${MAPTIME}\n" >> $MAPTIMESFILE
	 	 
	 # harvest toil stats
	 if test ${STATSFILE+defined}; then
		  echo "STATS FOR MAPEVAL ${JOB_STORE}-${REGION,,}-${GRAPH}" >> ${STATSFILE}
		  toil stats ${JOB_STORE}-${REGION,,} >> ${STATSFILE}
		  printf "\n\n" >> ${STATSFILE}
	 fi

	 toil clean ${JOB_STORE}-${REGION,,}
	 
	 # harvest the stats output and append it to our table
	 rm -f temp_stats.tsv
	 $CP_CMD $STATS_TSV ./temp_stats.tsv
	 if [ -f "./temp_stats.tsv" ] ; then
		  printf "${REGION}\t${GRAPH}\t$(cat ./temp_stats.tsv)\n" >> ${MEFILE}
	 else
		  printf "${REGION}\t${GRAPH}\tERROR\n" >> ${MEFILE}
	 fi	 
}	 

# grab the total time spent on vg map from the log
function extract_map_time {
	 local LOGFILE=$1

	 # we expect something along the lines of Successfully docker ran vg map graph.vg .... in XXX seconds
	 # so we extract and total all the XXXs.
	 grep "Successfully " $LOGFILE  | grep "ran vg map" | awk '{sum+= $(NF-1)}; END {print sum}'
}
	 
# Run the regions in series

rm -f $F1FILE
if test ${MEFILE+defined}; then
	 rm -rf $MEFILE
fi

run_bakeoff_region BRCA1 17 43044293 snp1kg
if [ "$FAST" != "1" ]; then
	 run_bakeoff_region BRCA2 13 32314860 snp1kg
	 run_bakeoff_region SMA 5 69216818 snp1kg
	 run_bakeoff_region LRC-KIR 19 54025633 snp1kg
	 run_bakeoff_region MHC 6 28510119 snp1kg
fi

echo "Created the following:"
echo ${JOB_STORE}-brca1-snp1kg
echo ${OUT_STORE}-brca1-snp1kg
if test ${MEFILE+defined}; then
	 echo ${OUT_STORE}-brca1-me-snp1kg
fi
if [ "$FAST" != "1" ]; then
	 echo ${JOB_STORE}-brca2-snp1kg
	 echo ${OUT_STORE}-brca2-snp1kg
	 echo ${JOB_STORE}-sma-snp1kg
	 echo ${OUT_STORE}-sma-snp1kg
	 echo ${JOB_STORE}-lrc-kir-snp1kg
	 echo ${OUT_STORE}-lrc-kir-snp1kg
	 echo ${JOB_STORE}-mhc-snp1kg
	 echo ${OUT_STORE}-mhc-snp1kg
	 if test ${MEFILE+defined}; then
		  echo ${OUT_STORE}-brca2-me-snp1kg
		  echo ${OUT_STORE}-sma-me-snp1kg
		  echo ${OUT_STORE}-lrc-kir-me-snp1kg
		  echo ${OUT_STORE}-mhc-me-snp1kg
	 fi
fi
echo "${F1FILE}"
        
