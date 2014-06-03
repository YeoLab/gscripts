__author__ = 'olga'

import unittest
from gscripts.mapping.repeat_align import RepeatAlign
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
        job_name = 'repeat_align'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        RepeatAlign(job_name=job_name, out_sh=out_sh,
                    directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N repeat_align
#PBS -o test_output/repeat_align.sh.out
#PBS -e test_output/repeat_align.sh.err
#PBS -V
#PBS -l walltime=2:30:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-6%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="bowtie         -c         -S         -q         -p 16         -e 100         -l 20         --un data/barcode_test.fastq.norep         all_ref         data/barcode_test.fastq         | grep -v "@"         |  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl         > data/barcode_test.fastq.repeat_counts"
cmd[2]="gunzip -c data/test.fastq.gz         |bowtie         -c         -S         -q         -p 16         -e 100         -l 20         --un data/test.fastq.gz.norep         all_ref         -         | grep -v "@"         |  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl         > data/test.fastq.gz.repeat_counts"
cmd[3]="gunzip -c data/test1_01_R1.fastq.gz         |bowtie         -c         -S         -q         -p 16         -e 100         -l 20         --un data/test1_01_R1.fastq.gz.norep         all_ref         -         | grep -v "@"         |  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl         > data/test1_01_R1.fastq.gz.repeat_counts"
cmd[4]="gunzip -c data/test1_01_R2.fastq.gz         |bowtie         -c         -S         -q         -p 16         -e 100         -l 20         --un data/test1_01_R2.fastq.gz.norep         all_ref         -         | grep -v "@"         |  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl         > data/test1_01_R2.fastq.gz.repeat_counts"
cmd[5]="gunzip -c data/test1_02_R1.fastq.gz         |bowtie         -c         -S         -q         -p 16         -e 100         -l 20         --un data/test1_02_R1.fastq.gz.norep         all_ref         -         | grep -v "@"         |  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl         > data/test1_02_R1.fastq.gz.repeat_counts"
cmd[6]="gunzip -c data/test1_02_R2.fastq.gz         |bowtie         -c         -S         -q         -p 16         -e 100         -l 20         --un data/test1_02_R2.fastq.gz.norep         all_ref         -         | grep -v "@"         |  perl /home/ppliu/tscc_scripts/count_aligned_from_sam.pl         > data/test1_02_R2.fastq.gz.repeat_counts"
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