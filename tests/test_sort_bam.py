__author__ = 'olga'

import unittest
from gscripts.mapping.sort_bam import SortBam
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

    def test_sort_bam(self):
        job_name = 'sort_bam'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        SortBam(job_name=job_name, out_sh=out_sh,
                directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N sort_bam
#PBS -o test_output/sort_bam.sh.out
#PBS -e test_output/sort_bam.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=8
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-4%10

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="samtools sort -@ 8 -m 50000000000 data/test.bam data/test.bam.sorted"
cmd[2]="samtools sort -@ 8 -m 50000000000 data/test.sorted.bam data/test.sorted.bam.sorted"
cmd[3]="samtools sort -@ 8 -m 50000000000 data/test_barcode_collapse.bam data/test_barcode_collapse.bam.sorted"
cmd[4]="samtools sort -@ 8 -m 50000000000 data/test_barcode_collapse.bc.bam data/test_barcode_collapse.bc.bam.sorted"
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