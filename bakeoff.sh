#!/usr/bin/env bash

## TODO: Get this going in some sort of CI framework
##       Ideally:  vg merge to master triggers docker build then a run of this script
##                 with the accuracy results logged to some kind of persistent leaderboard
##
## For now, all we have is the basic logic to do one run of bakeoff, expecting
## virtualenv etc has already been set up. 

if [ "$#" -ne 2 ]; then
	 echo "Syntax $0 <MESOS> <F1FILE>"
	 echo "(MESOS = [0,1])"
	 exit 1
fi

# all input data expected to be here:
BAKEOFF_STORE="s3://glennhickey-bakeoff-store"

# set to 1 to toggle mesos / s3
MESOS=$1
F1FILE=$2

# General Options
OPTS='--index_cores 7 --alignment_cores 7 --calling_cores 7 --vcfeval_cores 7 --realTimeLogging --logInfo'

# Hack in support for switching between mesos and local here
# (Note, for job store and out store, we will tack on -REGION to make them unique)
if [ "$MESOS" == "1" ]; then
	 JOB_STORE="aws:us-west-2:glennhickey-bakeoff-job-store"
	 OUT_STORE="aws:us-west-2:glennhickey-bakeoff-out-store"
	 CP_OUT_STORE="s3://glennhickey-bakeoff-out-store"
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
	 toil-vg run ${JOB_STORE}-${REGION,,} ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.fq.gz NA12878 ${OUT_STORE}-${REGION,,} --chroms ${CHROM} --call_opts "--offset ${OFFSET}" --graphs ${BAKEOFF_STORE}/snp1kg-${REGION}.vg --vcfeval_baseline ${BAKEOFF_STORE}/platinum_NA12878_${REGION}.vcf.gz --vcfeval_fasta ${BAKEOFF_STORE}/chrom.fa.gz ${BS_OPTS} ${OPTS}

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
run_bakeoff_region BRCA2 13 32314860
run_bakeoff_region SMA 5 69216818
run_bakeoff_region LRC_KIR 19 54025633
run_bakeoff_region MHC 6 28510119


        
