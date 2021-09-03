# Copyright (C) 2015 UCSC Computational Genomics Lab
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

define help

Supported targets: prepare, develop, sdist, clean, test, pypi.

Please note that all build targets require a virtualenv to be active.

The 'prepare' target installs toil-vg's build requirements into the current virtualenv.

The 'develop' target creates an editable install of toil-vg and its runtime requirements in the
current virtualenv. The install is called 'editable' because changes to the source code
immediately affect the virtualenv.

The 'sdist' target creates a source distribution of RNA-seq suitable for hot-deployment (not
implemented yet).

The 'clean' target undoes the effect of 'develop', and 'sdist'.

The 'test' target runs toil-vg's unit tests. Set the 'tests' variable to run a particular test, e.g.

	make test tests=src/toil_vg/test/test_vg.py::VGCGLTest::test_01_sim_small
	
By default it will try to use Docker. You can disable this by also passing container=None

	make test container=None tests=src/toil_vg/test/test_vg.py::VGCGLTest::test_01_sim_small

The 'pypi' target publishes the current commit of toil-vg to PyPI after enforcing that the working
copy and the index are clean, and tagging it as an unstable .dev build.

endef
export help
help:
	@echo "$$help"


python=python3
pip=pip3
tests=src
container=Docker
extras=


ifeq ($(shell uname -s),Darwin)
	# The Scipy build assumes it is working with Official Apple Clang as the default compiler.
	# But people building VG are likely to have GCC as the default compiler.
	# So make sure to build all Python modules with Clang, which they expect.
	export CC=/usr/bin/clang
	export CXX=/usr/bin/clang++
endif

green=\033[0;32m
normal=\033[0m
red=\033[0;31m


develop: check_venv
	$(pip) install -e .$(extras)
clean_develop: check_venv
	- $(pip) uninstall -y toil
	- rm -rf src/*.egg-info


sdist: check_venv
	$(python) setup.py sdist
clean_sdist:
	- rm -rf dist


test: check_venv check_build_reqs
	TOIL_VG_TEST_CONTAINER=$(container) $(python) setup.py test --pytest-args "-s -vv $(tests) --junitxml=test-report.xml"

pypi: check_venv check_clean_working_copy
	test "$$CI_COMMIT_REF_NAME" != "master" \
	&& echo "We're building a PR, skipping PyPI." || ( \
	set -x \
	&& tag_build=`$(python) -c 'pass;\
		from version import version as v;\
		from pkg_resources import parse_version as pv;\
		import os;\
		print("--tag-build=.dev" + os.getenv("CI_PIPELINE_IID") if pv(v).is_prerelease else "")'` \
	&& $(python) setup.py egg_info $$tag_build sdist bdist_egg upload )
clean_pypi:
	- rm -rf build/

clean: clean_develop clean_sdist clean_pypi clean_prepare


check_build_reqs:
	@$(python) -c 'import pytest' \
		|| ( echo "$(red)Build requirements are missing. Run 'make prepare' to install them.$(normal)" ; false )


prepare: check_venv
	$(pip) install numpy
	# TODO scikit-learn can't even begin to install unless numpy is already there, so numpy has to be first and by itself.
	# See https://github.com/scikit-learn/scikit-learn/issues/4164
	$(pip) install scipy scikit-learn
	$(pip) install pytest 'toil[aws,mesos]==4.1.0' biopython pyvcf
	pip list
clean_prepare: check_venv
	$(pip) uninstall -y pytest biopython numpy scipy scikit-learn pyvcf

check_venv:
	@$(python) -c 'import sys, os; sys.exit( int( 0 if "VIRTUAL_ENV" in os.environ else 1 ) )' \
		|| ( echo "$(red)A virtualenv must be active.$(normal)" ; false )


check_clean_working_copy:
	@echo "$(green)Checking if your working copy is clean ...$(normal)"
	@git diff --exit-code > /dev/null \
		|| ( echo "$(red)Your working copy looks dirty.$(normal)" ; false )
	@git diff --cached --exit-code > /dev/null \
		|| ( echo "$(red)Your index looks dirty.$(normal)" ; false )
	@test -z "$$(git ls-files --other --exclude-standard --directory)" \
		|| ( echo "$(red)You have are untracked files:$(normal)" \
			; git ls-files --other --exclude-standard --directory \
			; false )

clean_docker:
	-cd docker && make clean

test_docker:
	cd docker && make test

docker:
	cd docker && make

push_docker: docker
	cd docker && make push

.PHONY: help \
		prepare \
		develop clean_develop \
		sdist clean_sdist \
		test \
		pypi clean_pypi \
		clean \
		check_venv \
		check_clean_working_copy \
		docker \
		push_docker \
		clean_docker \
		test_docker \
		check_build_reqs
