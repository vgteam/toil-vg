import logging
import shlex
import shutil
import subprocess
import tempfile
import textwrap
import filecmp
from contextlib import closing
from unittest import TestCase
from urlparse import urlparse
from uuid import uuid4

import os
import posixpath
from bd2k.util.iterables import concat
from boto.s3.connection import S3Connection, Bucket
from boto.s3.key import Key

log = logging.getLogger(__name__)


class VGCGLTest(TestCase):
    """
    These tests *can* be parameterized with the following optional environment variables:

    TOIL_SCRIPTS_TEST_TOIL_OPTIONS - a space-separated list of additional command line arguments to pass to Toil via
    the script entry point. Default is the empty string.

    TOIL_SCRIPTS_TEST_JOBSTORE - the job store locator to use for the tests. The default is a file: locator pointing
    at a local temporary directory.

    TOIL_SCRIPTS_TEST_NUM_SAMPLES - the number of sample lines to generate in the input manifest
    """

    @classmethod
    def setUpClass(cls):
        super(VGCGLTest, cls).setUpClass()
        # FIXME: pull up into common base class
        logging.basicConfig(level=logging.INFO)

    def setUp(self):
        self.input_dir = '/home/mesosbox/testvg'
        self.output_dir = '/home/mesosbox/testvg'
	self.sample_reads = self.input_dir+ '/NA12877.brca1.brca2.bam.fq'
        self.test_vg_graph = self.input_dir+ '/BRCA1_BRCA2_may6.vg'
        self.jobStore = 'aws:us-west-2:cmarkello-s3-testvg-jobstore'
        self.base_command = concat('toil-vg', '--batchSystem', 'mesos', '--mesosMaster=mesos-master:5050',
                                   '--realTimeLogging', '--edge_max', '5', '--kmer_size',
                                   '16', '--num_fastq_chunks', '3', '--call_chunk_size', '10000', '--overwrite',
                                   '--index_mode', 'gcsa-mem', '--include_primary',
                                   self.jobStore)
        if not os.path.exists('/home/mesosbox/testvg/'):
            os.mkdir('/home/mesosbox/testvg/')
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/BRCA1_BRCA2_may6.vg', self.test_vg_graph])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.brca1.brca2.bam.fq', self.sample_reads])
        subprocess.check_call(['toil', 'clean', self.jobStore])

    def test_chr13_sampleNA12877(self):
        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
                                   self.output_dir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '13', '--path_size', '84989')
        self._assertOutput('13.vcf')
 
    def test_chr17_sampleNA12877(self):
        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
                                   self.output_dir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '17', '--path_size', '81189')
        self._assertOutput('17.vcf')
    
    def test_chr13_17_sampleNA12877(self):
        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
                                   self.output_dir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '13', '17', '--path_size', '84989', '81189')
        self._assertOutput('NA12877.vcf')
    
    def _run(self, *args):
        args = list(concat(*args))
        log.info('Running %r', args)
        subprocess.check_call(args)

    def _assertOutput(self, testFile):
        subprocess.check_call(['aws', 's3', 'cp', 's3://cmarkello-vgtoil-test--files/normal_'+testFile, '/home/mesosbox/testvg/normal_'+testFile])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cmarkello-hgvmdebugtest-output/'+testFile+'.gz', '/home/mesosbox/testvg/'+testFile+'.gz'])
        subprocess.check_call(['gzip', '-df', '/home/mesosbox/testvg/' + testFile + '.gz'])
        self.assertTrue(filecmp.cmp('/home/mesosbox/testvg/normal_' + testFile, '/home/mesosbox/testvg/' + testFile))

    def tearDown(self):
        shutil.rmtree(self.input_dir)
        subprocess.check_call(['toil', 'clean', self.jobStore])
