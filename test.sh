#!/usr/bin/env bash
# test.sh: run the tests in a virtualenv

set -e

rm -Rf .venv
virtualenv --system-site-packages .venv
. .venv/bin/activate
make test
