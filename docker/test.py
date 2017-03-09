#!/usr/bin/env python2.7
import subprocess
import sys
import unittest


class TestVGPipeline(unittest.TestCase):

    def test_docker_call(self):
        # print sys.argv
        tool = ['quay.io/ucsc_cgl/toilvg-cgl-pipeline:{}'.format(tag)]
        base = ['docker', 'run']
        # Check base call for help menu
        out = check_docker_output(command=base + tool, assert_fail=False)
        self.assertTrue('Toil vg DNA-seq pipeline' in out)

    def test_docker_generate_config_help(self):
        # print sys.argv
        tool = ['quay.io/ucsc_cgl/toilvg-cgl-pipeline:{}'.format(tag), 'generate-config', '--help']
        base = ['docker', 'run']
        # Check base call for help menu
        out = check_docker_output(command=base + tool, assert_fail=False)
        self.assertTrue('usage: vg_evaluation_pipeline generate-config [-h] [--whole_genome]' in out)
        
    def test_docker_index_help(self):
        # print sys.argv
        tool = ['quay.io/ucsc_cgl/toilvg-cgl-pipeline:{}'.format(tag), 'index', '--help']
        base = ['docker', 'run']
        # Check base call for help menu
        out = check_docker_output(command=base + tool, assert_fail=False)
        self.assertTrue('usage: vg_evaluation_pipeline index [-h]' in out)
    
    def test_docker_map_help(self):
        # print sys.argv
        tool = ['quay.io/ucsc_cgl/toilvg-cgl-pipeline:{}'.format(tag), 'map', '--help']
        base = ['docker', 'run']
        # Check base call for help menu
        out = check_docker_output(command=base + tool, assert_fail=False)
        self.assertTrue('usage: vg_evaluation_pipeline map [-h]' in out)

    def test_docker_call_help(self):
        # print sys.argv
        tool = ['quay.io/ucsc_cgl/toilvg-cgl-pipeline:{}'.format(tag), 'call', '--help']
        base = ['docker', 'run']
        # Check base call for help menu
        out = check_docker_output(command=base + tool, assert_fail=False)
        self.assertTrue('usage: vg_evaluation_pipeline call [-h]' in out)
    
    def test_docker_vcfeval_help(self):
        # print sys.argv
        tool = ['quay.io/ucsc_cgl/toilvg-cgl-pipeline:{}'.format(tag), 'vcfeval', '--help']
        base = ['docker', 'run']
        # Check base call for help menu
        out = check_docker_output(command=base + tool, assert_fail=False)
        self.assertTrue('usage: vg_evaluation_pipeline vcfeval [-h]' in out)

def check_docker_output(command, assert_fail=True):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = process.communicate()
    if assert_fail:
        assert process.returncode == 1
    else:
        assert process.returncode == 0
    return output[0]


if __name__ == '__main__':
    tag = sys.argv[1]
    del sys.argv[1]

    unittest.main()
