'''
Created on Aug 1, 2013

@author: Olga Botvinnik
'''
import unittest

from gscripts.qtools import Submitter
import subprocess
from subprocess import PIPE
import os
import shutil

import tests

HOSTNAME = subprocess.Popen('HOSTNAME', stdout=subprocess.PIPE).communicate()[
    0].strip()

# Global variable to test if we're on a server, e.g. TSCC or oolite (
# "compute" means one of the oolite compute nodes)
ON_SERVER = set([HOSTNAME]) & set(['tscc', 'oolite', 'compute'])


class Test(unittest.TestCase):
    commands = ['date', 'echo testing']
    out_dir = 'test_output'

    def setUp(self):
        os.mkdir(self.out_dir)

    def tearDown(self):
        shutil.rmtree(self.out_dir)

    def test_pbs(self):
        """Test PBS queue (TSCC)
        """
        job_name = 'test_qtools_submitter_pbs'
        submit_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        sub = Submitter(queue_type='PBS', sh_filename=submit_sh,
                        commands=self.commands,
                        job_name=job_name, nodes=1, ppn=1,
                        queue='home-yeo', walltime='0:01:00'
        )
        job_id = sub.job(submit=False)
        true_result_string = '''#!/bin/bash
#PBS -N test_qtools_submitter_pbs
#PBS -o {0}/test_qtools_submitter_pbs.sh.out
#PBS -e {0}/test_qtools_submitter_pbs.sh.err
#PBS -V
#PBS -l walltime=0:01:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home-yeo

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
date
echo testing
'''.format(self.out_dir)
        true_result = true_result_string.split('\n')

        # with open(submit_sh) as f:
        #     for x in f.readlines():
        #         print x,

        for true, test in zip(true_result, open(submit_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())

        # Make sure the job ID is a single (potentially multi-digit) integer
        # But only do this if we're on TSCC or oolite
        if ON_SERVER:
            self.assertRegexpMatches(job_id, '^\d+$')
            subprocess.Popen(["qdel", job_id],
                             stdout=PIPE)

    def test_sge(self):
        """Test SGE queue (oolite)
        """
        job_name = 'test_qtools_submitter_sge'
        submit_sh = '{}/{}.sh'.format(self.out_dir, job_name)
        sub = Submitter(queue_type='SGE', sh_filename=submit_sh,
                        commands=self.commands,
                        job_name=job_name, nodes=1, ppn=1,
                        queue='home-yeo', walltime='0:01:00'
        )
        job_id = sub.job(submit=False)
        true_result_string = '''#!/bin/bash
#$ -N test_qtools_submitter_sge
#$ -o {0}/test_qtools_submitter_sge.sh.out
#$ -e {0}/test_qtools_submitter_sge.sh.err
#$ -V
#$ -S /bin/bash
#$ -cwd
#$ -l bigmem
#$ -l h_vmem=16G
date
echo testing
'''.format(self.out_dir)
        true_result = true_result_string.split('\n')

        # with open(submit_sh) as f:
        #     for x in f.readlines():
        #         print x,
        for true, test in zip(true_result, open(submit_sh)):
            self.assertEqual(true.strip().split(), test.strip().split())

        # Make sure the job ID is a single (potentially multi-digit) integer
        # But only do this if we're on TSCC or oolite
        if ON_SERVER:
            self.assertRegexpMatches(job_id, '^\d+$')
            subprocess.Popen(["qdel", job_id],
                             stdout=PIPE)

#     def test_wait_for_pbs(self):
#         commands = ['date', 'echo testing PBS']
#         job_name = 'test_qtools_submitter_wait_for_pbs'
#         submit_sh = '%s.sh' % (job_name)
#         sub = Submitter(queue_type='PBS', sh_filename= submit_sh,
#                         commands=commands,
#                         job_name=job_name, wait_for=['11111'])
#         job_id = sub.write_sh(submit=False, nodes=1, ppn=16,
#                               queue='home-yeo', walltime='0:01:00')
#         true_result_string = '''#!/bin/bash
# #PBS -N test_qtools_submitter_wait_for_pbs
# #PBS -o test_qtools_submitter_wait_for_pbs.sh.out
# #PBS -e test_qtools_submitter_wait_for_pbs.sh.err
# #PBS -V
# #PBS -l walltime=0:01:00
# #PBS -l nodes=1:ppn=16
# #PBS -A yeo-group
# #PBS -q home-yeo
# #PBS -W depend=afterok:11111
#
# # Go to the directory from which the script was called
# cd $PBS_O_WORKDIR
# date
# echo testing PBS
# '''
#         true_result = true_result_string.split('\n')
#
#         # with open(submit_sh) as f:
#         #     for x in f.readlines():
#         #         print x,
#
#         for true, test in zip(true_result, open(submit_sh)):
#             self.assertEqual(true.strip().split(), test.strip().split())
#
#         # Make sure the job ID is a single (potentially multi-digit) integer
#         if ON_SERVER:
#             self.assertRegexpMatches(job_id, '^\d+$')
#             subprocess.Popen(["qdel", job_id], stdout=PIPE)
#
#     def test_wait_for_array_pbs(self):
#         commands = ['date', 'echo testing PBS']
#         job_name = 'test_qtools_submitter_wait_for_pbs'
#         submit_sh = '%s.sh' % (job_name)
#         sub = Submitter(queue_type='PBS', sh_filename= submit_sh,
#                         commands=commands,
#                         job_name=job_name, wait_for_array=['11111'])
#         job_id = sub.write_sh(submit=False, nodes=1, ppn=16,
#                               queue='home-yeo', walltime='0:01:00')
#         true_result_string = '''#!/bin/bash
# #PBS -N test_qtools_submitter_wait_for_pbs
# #PBS -o test_qtools_submitter_wait_for_pbs.sh.out
# #PBS -e test_qtools_submitter_wait_for_pbs.sh.err
# #PBS -V
# #PBS -l walltime=0:01:00
# #PBS -l nodes=1:ppn=16
# #PBS -A yeo-group
# #PBS -q home-yeo
# #PBS -W depend=afterokarray:11111
#
# # Go to the directory from which the script was called
# cd $PBS_O_WORKDIR
# date
# echo testing PBS
# '''
#         true_result = true_result_string.split('\n')
#
#         # with open(submit_sh) as f:
#         #     for x in f.readlines():
#         #         print x,
#
#         for true, test in zip(true_result, open(submit_sh)):
#             self.assertEqual(true.strip().split(), test.strip().split())
#
#         # Make sure the job ID is a single (potentially multi-digit) integer
#         if ON_SERVER:
#             self.assertRegexpMatches(job_id, '^\d+$')
#             subprocess.Popen(["qdel", job_id], stdout=PIPE)


if __name__ == "__main__":
    # import sys;sys.argv = ['', 'Test.test_main']
    unittest.main()