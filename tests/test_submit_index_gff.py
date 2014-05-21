__author__ = 'olga'

import unittest
from gscripts.miso.submit_index_gff import IndexGFF
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

    def test_submit_index_gff(self):
        job_name = 'index_gff'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        gff = '/projects/ps-yeolab/genomes/hg19/miso/SE.hg19.gff3'
        index_dir = '/projects/ps-yeolab/genomes/hg19/miso/SE_index'
        IndexGFF(gtf=None, gff=gff, index_dir=index_dir,
                 job_name=job_name, out_sh=out_sh,
                 submit=False)
        true_result = """#!/bin/bash
#PBS -N index_gff
#PBS -o test_output/index_gff.sh.out
#PBS -e test_output/index_gff.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
index_gff.py --index /projects/ps-yeolab/genomes/hg19/miso/SE.hg19.gff3 /projects/ps-yeolab/genomes/hg19/miso/SE_index
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())

    def test_submit_index_gff_gtf(self):
        job_name = 'index_gff'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        gtf = '/projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19' \
              '.annotation.gtf'
        index_dir = '/projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19' \
                    '.annotation_index'
        IndexGFF(gtf=gtf, gff=None, index_dir=index_dir,
                 job_name=job_name, out_sh=out_sh,
                 submit=False)
        true_result = """#!/bin/bash
#PBS -N index_gff
#PBS -o test_output/index_gff.sh.out
#PBS -e test_output/index_gff.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
gtf2gff3.pl /projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation.gtf > /projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation.gff
index_gff.py --index /projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation.gff /projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation_index
"""
        true_result = true_result.split('\n')
        # with open(out_sh) as f:
        #     for line in f:
        #         print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()