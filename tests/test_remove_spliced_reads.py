__author__ = 'olga'

import unittest
from gscripts.mapping.remove_spliced_reads import RemoveSpliced
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

    def test_removed_spliced(self):
        job_name = 'remove_spliced'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        RemoveSpliced(job_name=job_name, out_sh=out_sh,
                      directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N remove_spliced
#PBS -o test_output/remove_spliced.sh.out
#PBS -e test_output/remove_spliced.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-4%10

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="samtools view -h -F 4 data/test.bam | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > data/test.bam.unspliced.bam"
cmd[2]="samtools view -h -F 4 data/test.sorted.bam | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > data/test.sorted.bam.unspliced.bam"
cmd[3]="samtools view -h -F 4 data/test_barcode_collapse.bam | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > data/test_barcode_collapse.bam.unspliced.bam"
cmd[4]="samtools view -h -F 4 data/test_barcode_collapse.bc.bam | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > data/test_barcode_collapse.bc.bam.unspliced.bam"
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