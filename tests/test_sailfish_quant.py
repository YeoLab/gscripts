__author__ = 'olga'

import unittest
from gscripts.rnaseq.sailfish_quant import SailfishQuant
import tests
import os
import shutil
import sys


class Test(unittest.TestCase):
    out_dir = 'test_output'

    def setUp(self):
        os.mkdir(self.out_dir)

    def tearDown(self):
        shutil.rmtree(self.out_dir)

    def test_sailfish_quant_single_end(self):
        job_name = 'sailfish_quant'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        SailfishQuant('read1.fastq.gz', None, 'sailfish_index', 'sailfish_out',
                      job_name=job_name, out_sh=out_sh, submit=False)
        true_result = """#!/bin/bash
#PBS -N sailfish_quant
#PBS -o test_output/sailfish_quant.sh.out
#PBS -e test_output/sailfish_quant.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=8
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
sailfish quant --index sailfish_out -l T=SE:S=U -1 <(gunzip read1.fastq.gz)  --out sailfish_index --threads 8
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())

    def test_sailfish_quant_paired_end(self):
        job_name = 'sailfish_quant'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        SailfishQuant('read1.fastq.gz', 'read2.fastq.gz', 'sailfish_index',
                      'sailfish_out',
                      job_name=job_name, out_sh=out_sh, submit=False)
        true_result = """#!/bin/bash
#PBS -N sailfish_quant
#PBS -o test_output/sailfish_quant.sh.out
#PBS -e test_output/sailfish_quant.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=8
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
sailfish quant --index sailfish_out -l T=PE:O=><:S=SA -1 <(gunzip read1.fastq.gz) -2 <(gunzip read2.fastq.gz) --out sailfish_index --threads 8
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()