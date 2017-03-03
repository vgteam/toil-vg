#!/usr/bin/env bash

## TODO: Get this going in some sort of CI framework
##       Ideally:  vg merge to master triggers docker build then a run of this script
##                 with the accuracy results logged to some kind of persistent leaderboard
##
## For now, all we have is the basic logic to do one run of bakeoff, expecting
## virtualenv etc has already been set up. 

usage() { printf "Usage: $0 [Options] <Output-prefix> <Ouptut F1 File.tsv> \nOptions:\n\t-m\tmesos\n\t-f\tfast (just BRCA1)\n\t-D\tno Docker\n\t-t <N>\tthreads [DEFAULT=7]\n\t-c <F>\tconfig file\n\n" 1>&2; exit 1; }

FAST=0
MESOS=0
DOCKER=1
CORES=7
CONFIG=0

while getopts "fmDt:c:" o; do
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
OPTS="--gcsa_index_cores ${CORES} --kmers_cores ${CORES} --alignment_cores ${CORES} --calling_cores ${CORES} --vcfeval_cores ${CORES} --realTimeLogging --logInfo"
if [ "$CONFIG" != "0" ]; then
	 OPTS="${OPTS} --config ${CONFIG}"
fi

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
	 JOB_STORE="./bakeoff-job-store"
	 OUT_STORE="./bakeoff-out-store"
	 CP_OUT_STORE=${OUT_STORE}
	 BS_OPTS=""
	 CP_CMD="cp"
	 RM_CMD="rm"
fi

# Hack in support for turning Docker on and off
if [ "$DOCKER" == "0" ]; then
    DOCKER_OPTS="--no_docker"
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
	 toil-vg run ${JOB_STORE}-${REGION,,} NA12878 ${OUT_STORE}-${REGION,,} --fastq ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.fq.gz --chroms ${CHROM} --call_opts "--offset ${OFFSET}" --graphs ${BAKEOFF_STORE}/snp1kg-${REGION}.vg --vcfeval_baseline ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.vcf.gz --vcfeval_fasta ${BAKEOFF_STORE}/chrom.fa.gz ${BS_OPTS} ${DOCKER_OPTS} ${OPTS}

	 toil clean ${JOB_STORE}-${REGION,,}

	 # harvest the f1 output and append it to our table
	 rm -f temp_f1.txt
	 $CP_CMD $F1_SCORE ./temp_f1.txt
	 if [ -f "./temp_f1.txt" ] ; then
		  printf "${REGION}\t$(cat ./temp_f1.txt)\n" >> ${F1FILE}
	 else
		  printf "${REGION}\tERROR\n" >> ${F1FILE}
	 fi
}

# Run the regions in series

rm -f $F1FILE

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
if [ "$FAST" != "1" ]; then
	 echo ${JOB_STORE}-brca2
	 echo ${OUT_STORE}-brca2
	 echo ${JOB_STORE}-sma
	 echo ${OUT_STORE}-sma
	 echo ${JOB_STORE}-lrc-kir
	 echo ${OUT_STORE}-lrc-kir
	 echo ${JOB_STORE}-mhc
	 echo ${OUT_STORE}-mhc
fi

        
