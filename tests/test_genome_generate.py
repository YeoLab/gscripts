__author__ = 'olga'

import unittest
from gscripts.mapping.genome_generate import GenomeGenerate
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

    def test_genome_generate(self):
        true_result = """#!/bin/bash
#PBS -N genome_generate
#PBS -o test_output/genome_generate.sh.out
#PBS -e test_output/genome_generate.sh.err
#PBS -V
#PBS -l walltime=4:00:00
#PBS -l nodes=1:ppn=16
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
STAR --runMode genomeGenerate --genomeDir /projects/ps-yeolab/genomes/hg19/star_sjdb/ --genomeFastaFiles /projects/ps-yeolab/genomes/hg19/gencode/v19/GRCh37.p13.genome.fa --runThreadN 16 --sjdbGTFfile /projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation.gtf --sjdbOverhang 100
"""

        true_result = true_result.split('\n')

        out_sh = '{}/genome_generate.sh'.format(self.out_dir)
        genomeFastaFiles = '/projects/ps-yeolab/genomes/hg19/gencode/v19/GRCh37.p13.genome.fa'
        sjdb = '--sjdbGTFfile /projects/ps-yeolab/genomes/hg19/gencode/v19/gencode.v19.annotation.gtf'
        GenomeGenerate(genomeDir='/projects/ps-yeolab/genomes/hg19/star_sjdb/',
                       genomeFastaFiles=genomeFastaFiles,
                       sjdb=sjdb, sjdbOverhang=100,
                       job_name='genome_generate',
                       out_sh=out_sh,
                       submit=False)
        # with open(out_sh) as f:
        #     for line in f:
        #         sys.stdout.write(line)

        for true, test in zip(true_result, open(out_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())


if __name__ == "__main__":
    unittest.main()