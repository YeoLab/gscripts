__author__ = 'olga'

import unittest
from gscripts.rnaseq.rpkm_from_counts import RPKMfromCounts
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

    def test_rpkm_from_counts(self):
        job_name = 'rpkm_from_counts'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        RPKMfromCounts(job_name=job_name, out_sh=out_sh,
                       directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N rpkm_from_counts
#PBS -o test_output/rpkm_from_counts.sh.out
#PBS -e test_output/rpkm_from_counts.sh.err
#PBS -V
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-2%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="single_RPKM.py -i data/test.count -o data/test.count.rpkm"
cmd[2]="single_RPKM.py -i data/test_single_RPKM.count -o data/test_single_RPKM.count.rpkm"
eval ${cmd[$PBS_ARRAYID]}
"""
        true_result = true_result.split('\n')
        with open(out_sh) as f:
            for line in f:
                print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()