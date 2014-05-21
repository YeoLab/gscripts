__author__ = 'olga'

import unittest
from gscripts.rnaseq.sailfish_index import SailfishIndex
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

    def test_sailfish_index(self):
        job_name = 'sailfish_index'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        SailfishIndex(fasta='data/test.fasta', kmer_size=31,
                      job_name=job_name, out_sh=out_sh, submit=False)
        true_result = """#!/bin/bash
#PBS -N sailfish_index
#PBS -o test_output/sailfish_index.sh.out
#PBS -e test_output/sailfish_index.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=8
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
sailfish index --kmerSize 31 --threads 8 --transcripts data/test.fasta --out data/test.fasta_sailfish_index_k31
"""
        true_result = true_result.split('\n')
        with open(out_sh) as f:
            for line in f:
                print line,

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()