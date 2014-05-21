__author__ = 'olga'

import unittest
from gscripts.rnaseq.count_bam import CountBam
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

    def test_count_bam(self):
        job_name = 'count_bam'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        CountBam(orientation='both', bed='data/test.bed',
                 job_name=job_name, out_sh=out_sh,
                 directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N count_bam
#PBS -o test_output/count_bam.sh.out
#PBS -e test_output/count_bam.sh.err
#PBS -V
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-1%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="count_tags.py --annotation_file data/test.bed -f both -b data/test.sorted.bam -o data/test.sorted.bam.count"
eval ${cmd[$PBS_ARRAYID]}
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()