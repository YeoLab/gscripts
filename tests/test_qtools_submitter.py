'''
Created on Aug 1, 2013

@author: Olga Botvinnik
'''
import unittest

from qtools import Submitter

import tests


class Test(unittest.TestCase):
    def test_main(self):
        commands = ['date', 'echo testing PBS']
        job_name = 'test_qtools_submitter'
        sub = Submitter(queue_type='PBS', sh_file='%s.sh' % job_name,
                        command_list=commands,
                        job_name=job_name)
        job_id = sub.write_sh(submit=True, nodes=1, ppn=16,
                                 queue='home-yeo', walltime='0:01:00')

        for true, test in zip(true_result, open(tests.get_file(rpkm_file))):
            self.assertEqual(true.strip().split(), test.strip().split())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_main']
    unittest.main()