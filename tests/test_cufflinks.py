__author__ = 'olga'

import unittest
from gscripts.rnaseq.cufflinks import Cufflinks
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

    def test_cufflinks(self):
        job_name = 'cufflinks'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        Cufflinks(gtf='data/test.gtf',
                  job_name=job_name, out_sh=out_sh,
                  directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N cufflinks
#PBS -o test_output/cufflinks.sh.out
#PBS -e test_output/cufflinks.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=8
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-1%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="cufflinks --GTF data/test.gtf --GTF-guide  --multi-read-correct --num-threads 8 data/test.sorted.bam"
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