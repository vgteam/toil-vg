import sys

from version import version, required_versions
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
    tests_require=['pytest==2.8.3'],
    package_dir={'': 'src'},
    packages=find_packages('src'),
    entry_points={
        'console_scripts': ['toil-vg = toil_vg.vg_evaluation_pipeline:main',
                            'toil-vg-call = toil_vg.vg_call:main',
                            'toil-vg-index = toil_vg.vg_index:main',
                            'toil-vg-map = toil_vg.vg_map:main']})

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
      "pip install toil[aws,mesos]%s"
      "\n\n"
      "on every EC2 instance. Refer to Toil's documentation at http://toil.readthedocs.io/en/latest/installation.html "
      "for more information."
      % required_versions['toil'])
