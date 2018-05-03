#!/usr/bin/env bash
# convenience script to run toil-vg commands on a cluster created via create-ec2-leader.sh
# pass through toil-vg arguments but not "toil-vg" itself or any toil options. 

# Need Toil installed
if ! [ -x "$(command -v toil)" ]; then
	 printf "Toil must be installed in order to create leader nodes\n"
	 printf "ex:\n"
	 printf "virtualenv venv\n"
	 printf ". venv/bin/activate\n"
	 printf "pip install -U toil[aws,mesos]==3.11.0\n"
	 exit 1
fi

NODE_TYPE="r3.8xlarge:0.85"
MAX_NODES="8"
usage() {
    # Print usage to stderr
    exec 1>&2
    printf "Usage: $0 [Options] QUOTED-TOIL-VG-ARGS \n"
    printf "Options:\n"
	 printf "    -n NODE_TYPE use this node_type (default ${NODE_TYPE})\n"
	 printf "    -m MAX_NODES use this max_nodes (default ${MAX_NODES})\n"
	 printf "    -l CLUSTER_NAME run on this leader via toil ssh-cluster (instead of current machine)\n"
	 printf "\n"
    exit 1
}

while getopts "hn:m:l:" o; do
    case "${o}" in
        n)
            NODE_TYPE=$OPTARG
            ;;
        m)
            MAX_NODES=$OPTARG
            ;;
        l)
            CLUSTER_NAME=$OPTARG
            ;;
		  
        *)
            usage
            ;;
    esac
done

shift $((OPTIND-1))

if [[ "$#" -lt "1" ]]; then
    # Too few arguments
    usage
fi

TOIL_VG_ARGS="${1}"
shift

set -x

# Some default Toil options.  Need to edit in here to change (or add cli option)
# In particular, gcsa indexing of whole genome may need more disk than 3.8xlarge
NODE_OPTS="--nodeTypes ${NODE_TYPE} --maxNodes ${MAX_NODES} --defaultPreemptable"
RETRY_OPTS="--retryCount 3"
LOG_OPTS="--realTimeLogging --logInfo --realTimeStderr"
TOIL_VG_OPTS=""
# We need the master's IP to make Mesos go
if [ -z ${CLUSTER_NAME} ]
then
	 MASTER_IP="$(hostname -i)"
else
	 MASTER_IP="$($PREFIX toil ssh-cluster --insecure --zone=us-west-2a --logOff "${CLUSTER_NAME}" hostname -i)"
fi
MASTER_IP="$(printf $MASTER_IP | sed -e 's/\r$//' -e 's/\n//')"
MESOS_OPTS="--batchSystem=mesos --mesosMaster=${MASTER_IP}:5050"
# Put together our Toil Options
TOIL_OPTS="--provisioner aws ${NODE_OPTS} ${RETRY_OPTS} ${LOG_OPTS} ${MESOS_OPTS} ${TOIL_VG_OPTS}"

if [ -z ${CLUSTER_NAME} ]
then
	 EVAL_PREFIX="toil-vg"
else
	 EVAL_PREFIX="eval $PREFIX toil ssh-cluster --insecure --zone=us-west-2a "${CLUSTER_NAME}" venv/bin/toil-vg"
fi

# Run our toil command
${EVAL_PREFIX} ${TOIL_VG_ARGS} ${TOIL_OPTS}
TOIL_ERROR="$?"

