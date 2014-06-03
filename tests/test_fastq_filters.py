__author__ = 'olga'

import unittest
from gscripts.mapping.fastq_filters import FastqFilters
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

    def test_fastq_filters(self):
        true_result = """#!/bin/bash
#PBS -N fatsq_filters
#PBS -o test_output/fastq_filters.sh.out
#PBS -e test_output/fastq_filters.sh.err
#PBS -V
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=2
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-1%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="echo data/test.fastq.gz; zcat data/test.fastq.gz | fastx_artifacts_filter | fastq_quality_trimmer -l 20 -t 30 | fastq_quality_filter -q 30 -p 90 -z > filtered/data/test.fastq.gz"
eval ${cmd[$PBS_ARRAYID]}
"""
        true_result = true_result.split('\n')

        out_sh = '{}/fastq_filters.sh'.format(self.out_dir)
        FastqFilters(job_name='fatsq_filters', submit=False,
                     directory='data/', out_sh=out_sh)
        # with open(out_sh) as f:
        #     for line in f:
        #         sys.stdout.write(line)

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()