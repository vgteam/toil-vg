import logging
import shlex
import shutil
import subprocess
import tempfile
import textwrap
import filecmp
import vcf
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
        self.workdir = tempfile.mkdtemp()
        self.jobStoreAWS = 'aws:us-west-2:testvg-{}'.format(uuid4())
        self.jobStoreLocal = '{}/local-testvg-{}'.format(self.workdir, uuid4())
        self.base_command = concat('toil-vg',
                                   '--realTimeLogging', '--logDebug', '--edge_max', '5', '--kmer_size',
                                   '16', '--num_fastq_chunks', '4', '--call_chunk_size', '20000', '--overwrite',
                                   '--index_mode', 'gcsa-mem', '--include_primary', '--index_cores', '8', '--alignment_cores', '8',
                                   '--calling_cores', '8')
        
        # default output store
        self.outstore = 'aws:us-west-2:cmarkello-hgvmdebugtest-output'
        self.local_outstore = os.path.join(self.workdir, 'cmarkello-hgvmdebugtest-output')

        # copy the sequence information for vcf comparison
        # (lumped in one file out of laziness.  todo: at least split by chromosome)
        self.chrom_fa = os.path.join(self.workdir, 'chrom.fa.gz')
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/chrom.fa.gz', self.chrom_fa])
        
    def test_chr17_sampleNA12877(self):
        ''' Test sample BRCA1 output, graph construction and use, and local file processing
        '''
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.brca1.bam.fq', self.workdir])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/BRCA1_chrom_name_chop_100.vg', self.workdir])
        self.sample_reads = os.path.join(self.workdir, 'NA12877.brca1.bam.fq')
        self.test_vg_graph = os.path.join(self.workdir, 'BRCA1_chrom_name_chop_100.vg')
        self._run(self.base_command, self.jobStoreLocal, self.test_vg_graph, self.sample_reads, 'NA12877', self.workdir,
                                   self.local_outstore, '--path_name', '17', '--path_size', '81189', '--call_opts', '--offset 43044293' )

        self._assertOutput('17.vcf', self.local_outstore)

    def test_chr19_sampleNA12877(self):
        ''' Test sample LRC KIR output
        '''
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.lrc_kir.bam.small.fq', self.workdir])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/LRC_KIR_chrom_name_chop_100.small.vg', self.workdir])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/lrc_kir_index.tar.gz', self.workdir])
        self.sample_reads = os.path.join(self.workdir, 'NA12877.lrc_kir.bam.small.fq')
        self.test_vg_graph = os.path.join(self.workdir, 'LRC_KIR_chrom_name_chop_100.small.vg')
        self.test_index = os.path.join(self.workdir, 'lrc_kir_index.tar.gz')
        self._run(self.base_command, self.jobStoreAWS, self.test_vg_graph, self.sample_reads, 'NA12877', self.workdir,
                                   self.outstore,  '--gcsa_index', self.test_index, '--path_name', '19', '--path_size', '50000')
        self._assertOutput('19.vcf', self.outstore)

    def test_chr6_MHC_sampleNA12877(self):
        ''' Test sample MHC output
        '''
        self.sample_reads = 's3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.mhc.bam.small.fq'
        self.test_vg_graph = 's3://cgl-pipeline-inputs/vg_cgl/ci/MHC_chrom_name_chop_100.small.vg'
        self.test_index = 's3://cgl-pipeline-inputs/vg_cgl/ci/mhc_index.tar.gz'
        self._run(self.base_command, self.jobStoreAWS, self.test_vg_graph, self.sample_reads, 'NA12877', self.workdir,
                                   self.outstore,  '--gcsa_index', self.test_index, '--path_name', '6', '--path_size', '50000')
        self._assertOutput('6.vcf', self.outstore)

    def test_chr5_SMA_sampleNA12877(self):
        ''' Test sample SMA output
        '''
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/NA12877.sma.bam.small.fq', self.workdir])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/SMA_chrom_name_chop_100.small.vg', self.workdir])
        subprocess.check_call(['aws', 's3', 'cp', 's3://cgl-pipeline-inputs/vg_cgl/ci/sma_index.tar.gz', self.workdir])
        self.sample_reads = os.path.join(self.workdir, 'NA12877.sma.bam.small.fq')
        self.test_vg_graph = os.path.join(self.workdir, 'SMA_chrom_name_chop_100.small.vg')
        self.test_index = os.path.join(self.workdir, 'sma_index.tar.gz')
        self._run(self.base_command, self.jobStoreAWS, self.test_vg_graph, self.sample_reads, 'NA12877', self.workdir,
                                   self.outstore, '--path_name', '5', '--path_size', '50000')
        self._assertOutput('5.vcf', self.outstore)
    
    def _run(self, *args):
        args = list(concat(*args))
        log.info('Running %r', args)
        subprocess.check_call(args)

    def _assertOutput(self, testFile, outstore, f1_threshold=0.90):
        s3_input_path = 's3://cgl-pipeline-inputs/vg_cgl/ci/'
        if outstore.startswith('aws:'):
            cp = ['aws', 's3', 'cp']
            output_path = outstore.replace('aws:us-west-2:','s3://')
        else:
            cp = ['cp']
            output_path = outstore
        testName = testFile + '.gz'
        normalName = 'normal_' + testFile + '.gz'
        localTest = os.path.join(self.workdir, testName)
        localNormal = os.path.join(self.workdir, normalName)

        # copy baseline
        subprocess.check_call(['aws', 's3', 'cp', os.path.join(s3_input_path, normalName), localNormal])
        subprocess.check_call(['aws', 's3', 'cp', os.path.join(s3_input_path, normalName) + '.tbi', localNormal + '.tbi'])

        # copy calls
        subprocess.check_call(cp + [os.path.join(output_path, testName), localTest])
        subprocess.check_call(cp + [os.path.join(output_path, testName) + '.tbi', localTest + '.tbi'])

        # run vcfeval
        eval_results = os.path.join(self.workdir, 'eval_results')
        subprocess.check_call(['toil-vg-vcfeval', self.jobStoreLocal, localTest, localNormal, self.chrom_fa, eval_results,
                               '--vcfeval_cores', '8'])
        with open(os.path.join(eval_results, 'f1.txt')) as f1_file:
            f1_score = float(f1_file.readline().strip())
        
        self.assertTrue(f1_score >= f1_threshold)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        subprocess.check_call(['toil', 'clean', self.jobStoreAWS])
        subprocess.check_call(['toil', 'clean', self.jobStoreLocal])
