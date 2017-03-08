#!/usr/bin/env bash



# Create Toil venv
rm -rf .env
virtualenv --never-download .env
. .env/bin/activate

# Prepare directory for temp files
TMPDIR=/mnt/ephemeral/tmp
rm -rf $TMPDIR
mkdir $TMPDIR
export TMPDIR

# Create s3am venv
rm -rf s3am
virtualenv --never-download s3am && s3am/bin/pip install s3am==2.0
mkdir -p bin
# Expose binaries to the PATH
ln -snf ${PWD}/s3am/bin/s3am bin/
export PATH=$PATH:${PWD}/bin

# Create awscli venv
rm -rf awscli
virtualenv --never-download awscli && awscli/bin/pip install awscli
# Expose binaries to the PATH
ln -snf ${PWD}/awscli/bin/aws bin/
export PATH=$PATH:${PWD}/bin

make prepare
make develop
make test
make docker
make test_docker
make clean

# clean working copy to satisfy corresponding check in Makefile
rm -rf bin awscli s3am
make pypi
make push_docker
make clean_docker

rm -rf .env $TMPDIR
