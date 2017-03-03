import logging
import shlex
import shutil
import subprocess
import tempfile
import textwrap
import filecmp
import vcf
from contextlib import closing
from unittest import TestCase, skip
from urlparse import urlparse
from uuid import uuid4

import os
import posixpath
from bd2k.util.iterables import concat
from boto.s3.connection import S3Connection, Bucket
from boto.s3.key import Key
from toil_vg.iostore import IOStore

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

    def _download_input(self, filename, local_path = None):
        tgt = os.path.join(self.workdir, filename) if local_path is None else local_path
        with open(tgt, 'w') as f:
            self.bucket.get_key('/vg_cgl/ci/{}'.format(filename)).get_contents_to_file(f)
        
    def setUp(self):
        self.workdir = tempfile.mkdtemp()
        self.jobStoreAWS = 'aws:us-west-2:testvg-{}'.format(uuid4())
        self.jobStoreLocal = '{}/local-testvg-{}'.format(self.workdir, uuid4())
        self.connection = S3Connection()
        self.bucket = self.connection.get_bucket('cgl-pipeline-inputs')
        self.base_command = concat('toil-vg', 'run',
                                   '--realTimeLogging', '--logInfo', '--reads_per_chunk', '8000',
                                   '--call_chunk_size', '20000',
                                   '--index_mode', 'gcsa-mem', '--gcsa_index_cores', '7', '--kmers_cores', '7',
                                   '--alignment_cores', '4',
                                   '--calling_cores', '4', '--vcfeval_cores', '4')
        
        # default output store
        self.outstore = 'aws:us-west-2:toilvg-jenkinstest-outstore-{}'.format(uuid4())
        self.local_outstore = os.path.join(self.workdir, 'toilvg-jenkinstest-outstore-{}'.format(uuid4()))

    def test_1_sim_small(self):
        ''' 
        This test uses simulated reads from the small dataset from vg, created as follows:
        vg construct -r test/small/x.fa -v test/small/x.vcf.gz > small.vg
        vg index -x small.xg small.vg
        vg sim -x small.xg  -l 100 -n 10000 -s 0 -e 0.001 -i 0.0001 > small_sim_reads
        # Reads converted to small_sim_reads.fq with script setting qualities to B
        (small.vcf.gz and small.fa.gz below are just x.vcf.gz and x.fa from input)
        '''
        self.sample_reads = 's3://cgl-pipeline-inputs/vg_cgl/ci/small_sim_reads.fq.gz'
        self.test_vg_graph = 's3://cgl-pipeline-inputs/vg_cgl/ci/small.vg'
        self.baseline = 's3://cgl-pipeline-inputs/vg_cgl/ci/small.vcf.gz'
        self.chrom_fa = 's3://cgl-pipeline-inputs/vg_cgl/ci/small.fa.gz'

        self._run(self.base_command, self.jobStoreLocal, 'sample',
                  self.local_outstore,  '--fastq', self.sample_reads,
                  '--graphs', self.test_vg_graph,
                  '--chroms', 'x', '--vcfeval_baseline', self.baseline,
                  '--vcfeval_fasta', self.chrom_fa, '--vcfeval_opts', ' --squash-ploidy')
        
        self._assertOutput('sample', self.local_outstore, f1_threshold=0.95)

    def test_2_sim_small_standalone(self):
        ''' 
        Same as above, but chain standalone tools instead of toil-vg run
        '''
        self.sample_reads = 's3://cgl-pipeline-inputs/vg_cgl/ci/small_sim_reads.fq.gz'
        self.test_vg_graph = 's3://cgl-pipeline-inputs/vg_cgl/ci/small.vg'
        self.baseline = 's3://cgl-pipeline-inputs/vg_cgl/ci/small.vcf.gz'
        self.chrom_fa = 's3://cgl-pipeline-inputs/vg_cgl/ci/small.fa.gz'
        
        self._run('toil-vg', 'index', self.jobStoreLocal, self.local_outstore,
                  '--graphs', self.test_vg_graph, '--chroms', 'x',
                   '--gcsa_index_cores', '7', '--kmers_cores', '7',
                  '--realTimeLogging', '--logInfo', '--index_name', 'small')

        self._run('toil-vg', 'map', self.jobStoreLocal, 'sample',
                  os.path.join(self.local_outstore, 'small.xg'),
                  os.path.join(self.local_outstore, 'small.gcsa'),
                  self.local_outstore,  '--fastq', self.sample_reads,
                  '--id_ranges', os.path.join(self.local_outstore, 'small_id_ranges.tsv'),
                  '--alignment_cores', '7', '--reads_per_chunk', '8000',
                  '--realTimeLogging', '--logInfo')
        
        self._run('toil-vg', 'call', self.jobStoreLocal,
                  os.path.join(self.local_outstore, 'small.xg'), 'sample',
                  self.local_outstore, '--gams', os.path.join(self.local_outstore, 'sample_x.gam'), 
                  '--chroms', 'x', '--call_chunk_size', '20000', '--calling_cores', '4',
                  '--realTimeLogging', '--logInfo')

        self._run('toil-vg', 'vcfeval', self.jobStoreLocal,
                  os.path.join(self.local_outstore, 'sample.vcf.gz'), self.baseline,
                  self.chrom_fa, self.local_outstore,
                  '--vcfeval_opts', ' --squash-ploidy',
                  '--realTimeLogging', '--logInfo')
        
        self._assertOutput(None, self.local_outstore, f1_threshold=0.95)

    def test_3_BRCA1_NA12877(self):
        ''' Test sample BRCA1 output, graph construction and use, and local file processing
        '''
        self._download_input('NA12877.brca1.bam_1.fq.gz')
        self._download_input('NA12877.brca1.bam_2.fq.gz')
        self._download_input('snp1kg-brca1.vg')
        self._download_input('platinum_NA12877_BRCA1.vcf.gz.tbi')
        self._download_input('platinum_NA12877_BRCA1.vcf.gz')
        self._download_input('BRCA1.fa.gz')
        
        self.sample_reads = os.path.join(self.workdir, 'NA12877.brca1.bam_1.fq.gz')
        self.sample_reads2 = os.path.join(self.workdir, 'NA12877.brca1.bam_2.fq.gz')
        self.test_vg_graph = os.path.join(self.workdir, 'snp1kg-brca1.vg')
        self.baseline = os.path.join(self.workdir, 'platinum_NA12877_BRCA1.vcf.gz')
        self.chrom_fa = os.path.join(self.workdir, 'BRCA1.fa.gz')
        
        self._run(self.base_command, self.jobStoreLocal, 'NA12877',
                  self.local_outstore, '--fastq', self.sample_reads, self.sample_reads2, '--graphs',
                  self.test_vg_graph, '--chroms', '17',
                  '--call_opts', '--offset 43044293',
                  '--vcfeval_baseline', self.baseline, '--vcfeval_fasta', self.chrom_fa)

        self._assertOutput('NA12877', self.local_outstore, f1_threshold=0.48)
    
    def _run(self, *args):
        args = list(concat(*args))
        log.info('Running %r', args)
        subprocess.check_call(args)

    def _assertOutput(self, sample_name, outstore, f1_threshold=0.90):

        # grab the f1.txt file
        local_f1 = os.path.join(self.workdir, 'f1.txt')
        io_store = IOStore.get(outstore)
        f1_path = 'vcfeval_output_f1.txt'
        if sample_name:
            f1_path = sample_name + '_' + f1_path
        io_store.read_input_file(f1_path, local_f1)

        with open(local_f1) as f1_file:
            f1_score = float(f1_file.readline().strip())
        self.assertTrue(f1_score >= f1_threshold)

    def tearDown(self):
        shutil.rmtree(self.workdir)
        subprocess.check_call(['toil', 'clean', self.jobStoreAWS])
        subprocess.check_call(['toil', 'clean', self.jobStoreLocal])
        
