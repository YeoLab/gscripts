__author__ = 'olga'

import unittest
from gscripts.rnaseq.submit_oldsplice import OldspliceSubmitter
import tests
import os
import shutil


class Test(unittest.TestCase):
    out_dir = 'test_output'

    def setUp(self):
        os.mkdir(self.out_dir)

    def tearDown(self):
        shutil.rmtree(self.out_dir)

    def test_submit_oldsplice(self):
        true_result = """#!/bin/bash
#PBS -N oldsplice
#PBS -o test_output/runOldsplice.sh.out
#PBS -e test_output/runOldsplice.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-4%1000

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="oldsplice.py -b test1.bam -s hg19 -o test1.splices --splice_type SE --splice_type MXE --processors 16"
cmd[2]="oldsplice.py -f -b test1.bam -s hg19 -o test1.flip.splices --splice_type SE --splice_type MXE --processors 16"
cmd[3]="oldsplice.py -b test2.bam -s hg19 -o test2.splices --splice_type SE --splice_type MXE --processors 16"
cmd[4]="oldsplice.py -f -b test2.bam -s hg19 -o test2.flip.splices --splice_type SE --splice_type MXE --processors 16"
eval ${cmd[$PBS_ARRAYID]}
"""
        true_result = true_result.split('\n')

        out_sh = '{}/runOldsplice.sh'.format(self.out_dir)
        OldspliceSubmitter(tests.get_file('sample_info.txt'), submit=False,
                           queue_type='PBS', out_sh=out_sh)
        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()