__author__ = 'olga'

import unittest
from gscripts.mapping.convert_sam import ConvertSam
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

    def test_convert_sam(self):
        job_name = 'convert_sam'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        ConvertSam(job_name=job_name, out_sh=out_sh,
                   queue_type='PBS', directory='data/',
                   submit=False)
        true_result = """#!/bin/bash
#PBS -N convert_sam
#PBS -o test_output/convert_sam.sh.out
#PBS -e test_output/convert_sam.sh.err
#PBS -V
#PBS -l walltime=1:00:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-1%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="samtools view -bS -q 10 data/test.sam > data/test.sam.bam"
eval ${cmd[$PBS_ARRAYID]}
"""
        true_result = true_result.split('\n')

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()