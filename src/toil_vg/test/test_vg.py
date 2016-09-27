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
#        self.workdir = '/home/mesosbox/testvg/'
        self.workdir = tempfile.mkdtemp()
        self.sample_reads = os.path.join(self.workdir, 'NA12877.brca1.brca2.bam.fq')
        self.test_vg_graph = os.path.join(self.workdir, 'BRCA1_BRCA2_may6.vg')
#        self.jobStore = '/var/lib/toil/cmarkello-s3-testvg-jobstore'
        self.jobStore = 'aws:us-west-2:testvg-{}'.format(uuid4())
#        self.jobStore = 'aws:us-west-2:cmarkello-s3-testvg-jobstore'
        self.base_command = concat('toil-vg',
                                   '--realTimeLogging', '--logDebug', '--edge_max', '5', '--kmer_size',
                                   '16', '--num_fastq_chunks', '3', '--call_chunk_size', '10000', '--overwrite',
                                   '--index_mode', 'gcsa-mem', '--include_primary', '--index_cores', '4', '--alignment_cores', '4',
                                   '--calling_cores', '4', self.jobStore)
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/BRCA1_BRCA2_may6.vg', self.test_vg_graph])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.brca1.brca2.bam.fq', self.sample_reads])
    
    def test_chr13_sampleNA12877(self):
        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
                                   self.workdir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '13', '--path_size', '84989')
        self._assertOutput('13.vcf')

#    def test_chr13_sampleNA12877(self):
#        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
#                                   self.workdir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
#                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '13', '--path_size', '84989')
#        self._assertOutput('13.vcf')
 
#    def test_chr17_sampleNA12877(self):
#        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
#                                   self.workdir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
#                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '17', '--path_size', '81189')
#        self._assertOutput('17.vcf')
    
#    def test_chr13_17_sampleNA12877(self):
#        self._run(self.base_command, self.test_vg_graph, self.sample_reads, 'NA12877',
#                                   self.workdir, 'aws:us-west-2:cmarkello-hgvmdebugtest-output',
#                                   'aws:us-west-2:cmarkello-hgvmdebugtest-output', '--path_name', '13', '17', '--path_size', '84989', '81189')
#        self._assertOutput('NA12877.vcf')
    
    def _run(self, *args):
        args = list(concat(*args))
        log.info('Running %r', args)
        subprocess.check_call(args)

    def _assertOutput(self, testFile):
        subprocess.check_call(['aws', 's3', 'cp', 's3://cmarkello-vgtoil-test--files/normal_'+testFile, self.workdir + 'normal_'+testFile])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cmarkello-hgvmdebugtest-output/'+testFile+'.gz', self.workdir + testFile+'.gz'])
        subprocess.check_call(['gzip', '-df', self.workdir + testFile + '.gz'])
        self.assertTrue(filecmp.cmp(self.workdir + 'normal_' + testFile, self.workdir +testFile))
#        subprocess.check_call(['aws', 's3', 'cp', 's3://cmarkello-vgtoil-test--files/normal_'+testFile, os.path.join(self.workdir, 'normal_'+testFile)])
#        subprocess.check_call(['aws', 's3', 'cp', 's3://cmarkello-hgvmdebugtest-output/'+testFile+'.gz', os.path.join(self.workdir, testFile+'.gz')])
#        subprocess.check_call(['gzip', '-df', os.path.join(self.workdir, testFile + '.gz')])
#        self.assertTrue(filecmp.cmp(os.path.join(self.workdir, 'normal_' + testFile), os.path.join(self.workdir, testFile)))

    def tearDown(self):
        shutil.rmtree(self.workdir)
        subprocess.check_call(['toil', 'clean', self.jobStore])
