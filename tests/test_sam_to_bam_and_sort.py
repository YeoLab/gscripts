__author__ = 'olga'

import unittest
from gscripts.mapping.sam_to_bam_and_sort import SamToBamAndSort
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

    def test_sam_to_bam_and_sort(self):
        job_name = 'sam_to_bam_and_sort'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        SamToBamAndSort(job_name=job_name, out_sh=out_sh,
                        directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N sam_to_bam_and_sort
#PBS -o test_output/sam_to_bam_and_sort.sh.out
#PBS -e test_output/sam_to_bam_and_sort.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
samtools view -bS data/test.sam > data/test.bam
samtools sort  data/test.bam data/test.sorted
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()