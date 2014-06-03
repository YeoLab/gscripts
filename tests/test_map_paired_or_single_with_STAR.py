__author__ = 'olga'

import unittest
from gscripts.mapping.map_paired_or_single_with_STAR import MapSTAR
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

    def test_map_paired_or_single_with_STAR(self):
        job_name = 'STAR'
        out_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        MapSTAR(genome='/projects/ps-yeolab/genomes/hg19/star_sjdb/',
                job_name=job_name, out_sh=out_sh,
                directory='data/', submit=False)
        true_result = """#!/bin/bash
#PBS -N STAR
#PBS -o test_output/STAR.sh.out
#PBS -e test_output/STAR.sh.err
#PBS -V
#PBS -l walltime=0:30:00
#PBS -l nodes=1:ppn=8
#PBS -A yeo-group
#PBS -q home
#PBS -t 1-2%20

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
cmd[1]="STAR         --runMode alignReads         --runThreadN 8         --genomeDir /projects/ps-yeolab/genomes/hg19/star_sjdb/         --genomeLoad LoadAndRemove         --readFilesCommand zcat         --readFilesIn           --outFileNamePrefix ./test1_01.         --outReadsUnmapped Fastx         --outFilterMismatchNmax 5         --outFilterMismatchNoverLmax 0.3         --outFilterMultimapNmax 5         --outFilterScoreMin 10         --outFilterType BySJout         --outSAMattributes All         --outSAMstrandField intronMotif         --clip5pNbases 0         --clip3pNbases 0         "
cmd[2]="STAR         --runMode alignReads         --runThreadN 8         --genomeDir /projects/ps-yeolab/genomes/hg19/star_sjdb/         --genomeLoad LoadAndRemove         --readFilesCommand zcat         --readFilesIn           --outFileNamePrefix ./test1_02.         --outReadsUnmapped Fastx         --outFilterMismatchNmax 5         --outFilterMismatchNoverLmax 0.3         --outFilterMultimapNmax 5         --outFilterScoreMin 10         --outFilterType BySJout         --outSAMattributes All         --outSAMstrandField intronMotif         --clip5pNbases 0         --clip3pNbases 0         "
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