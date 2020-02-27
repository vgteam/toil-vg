#!/usr/bin/env bash



# Create Toil venv
rm -rf .env
virtualenv -p python3.6 --never-download .env
. .env/bin/activate

# Upgrade pip3
pip3 install --upgrade pip3 setuptools

# Prepare directory for temp files
# Sometimes the instances have un-deletable files in tmp, so we continue through errors
TMPDIR=/mnt/ephemeral/tmp
set +e
rm -rf $TMPDIR
mkdir $TMPDIR
set -e
export TMPDIR

# Create s3am venv
rm -rf s3am
virtualenv -p python3.6 --never-download s3am && s3am/bin/pip3 install s3am==2.0
mkdir -p bin
# Expose binaries to the PATH
ln -snf ${PWD}/s3am/bin/s3am bin/
export PATH=$PATH:${PWD}/bin

# Create awscli venv
rm -rf awscli
virtualenv -p python3.6 --never-download awscli && awscli/bin/pip3 install awscli
# Expose binaries to the PATH
ln -snf ${PWD}/awscli/bin/aws bin/
export PATH=$PATH:${PWD}/bin

make prepare
make develop
make test
#make docker
#make test_docker
make clean

# clean working copy to satisfy corresponding check in Makefile
rm -rf bin awscli s3am
make pypi
#make push_docker
#make clean_docker

rm -rf .env $TMPDIR
