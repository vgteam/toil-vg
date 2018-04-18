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

	make test tests=src/toil/test/sort/sortTest.py::SortTest::testSort

The 'pypi' target publishes the current commit of toil-vg to PyPI after enforcing that the working
copy and the index are clean, and tagging it as an unstable .dev build.

endef
export help
help:
	@echo "$$help"


python=python2.7
pip=pip2.7
tests=src
extras=

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
	$(python) setup.py test --pytest-args "-vv $(tests) --junitxml=test-report.xml"

integration-test: check_venv check_build_reqs sdist
	TOIL_TEST_INTEGRATIVE=True $(python) run_tests.py integration-test $(tests)


pypi: check_venv check_clean_working_copy check_running_on_jenkins
	test "$$ghprbActualCommit" \
	&& echo "We're building a PR, skipping PyPI." || ( \
	set -x \
	&& tag_build=`$(python) -c 'pass;\
		from version import version as v;\
		from pkg_resources import parse_version as pv;\
		import os;\
		print "--tag-build=.dev" + os.getenv("BUILD_NUMBER") if pv(v).is_prerelease else ""'` \
	&& $(python) setup.py egg_info $$tag_build sdist bdist_egg upload )
clean_pypi:
	- rm -rf build/


clean: clean_develop clean_sdist clean_pypi clean_prepare


check_build_reqs:
	@$(python) -c 'import pytest' \
		|| ( echo "$(red)Build requirements are missing. Run 'make prepare' to install them.$(normal)" ; false )


prepare: check_venv
	$(pip) install numpy scipy scikit-learn==0.18.2
	$(pip) install pytest==2.8.3 'git+https://github.com/BD2KGenomics/toil.git#egg=toil[aws,mesos,azure]' biopython==1.67 pyvcf==0.6.8
	pip list
clean_prepare: check_venv
	$(pip) uninstall -y pytest biopython numpy scipy scikit-learn pyvcf

check_venv:
	@$(python) -c 'import sys; sys.exit( int( not hasattr(sys, "real_prefix") ) )' \
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


check_running_on_jenkins:
	@echo "$(green)Checking if running on Jenkins ...$(normal)"
	@test -n "$$BUILD_NUMBER" \
		|| ( echo "$(red)This target should only be invoked on Jenkins.$(normal)" ; false )

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
		check_running_on_jenkins \
		docker \
		push_docker \
		clean_docker \
		test_docker \
		check_build_reqs
