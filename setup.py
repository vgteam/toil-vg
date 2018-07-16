import sys
import os

# We can't import version.py because only one "version" module can ever be
# loaded in a Python process, and multiple setup.py scripts may have to run in
# the same process.
execfile(os.path.join(os.path.dirname(os.path.realpath(__file__)), "version.py"))

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand

kwargs = dict(
    name='toil-vg',
    version=version,
    description="UC Santa Cruz Computational Genomics Lab's Toil-based VG pipeline",
    author='UCSC Computational Genomics Lab',
    author_email='cgl-toil@googlegroups.com',
    url="https://github.com/BD2KGenomics/toil-vg",
    install_requires=[x + y for x, y in required_versions.iteritems()],
    dependency_links=[],
    tests_require=['pytest==2.8.3', 'numpy', 'scipy'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': ['toil-vg = toil_vg.vg_toil:main']}
)


class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        # Sanitize command line arguments to avoid confusing Toil code attempting to parse them
        sys.argv[1:] = []
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)

kwargs['cmdclass'] = {'test': PyTest}

setup(**kwargs)


print("\n\n"
      "Thank you for installing the UC Santa Cruz Computuational Genomics Lab VG DNA-seq pipeline! "
      "If you want to run this Toil-based pipeline on a cluster in a cloud, please install Toil "
      "with the appropriate extras. For example, To install AWS/EC2 support for example, run "
      "\n\n"
      "pip install \'toil[aws,mesos]>=3.6.0\'"
      "\n\n"
      "on every EC2 instance. Refer to Toil's documentation at http://toil.readthedocs.io/en/latest/installation.html "
      "for more information.")
