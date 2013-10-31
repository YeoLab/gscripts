'''
Created on Jun 21, 2013

@author: gabrielp
'''
import os
import unittest

from gscripts.rnaseq import single_RPKM

import tests


class Test(unittest.TestCase):

    def test_main(self):
        rpkm_file = "rpkm_test.rpkm"
        single_RPKM.main(tests.get_file("test_single_RPKM.count"), os.path.join(tests.get_test_dir(), rpkm_file ))
        
        true_result = ["gene    flag    RPKM",
                        "ENSG1    0    5025125.62814",
                         "ENSG2    0    0.0"]
                         
        for true, test in zip(true_result, open(tests.get_file(rpkm_file))):
            self.assertEqual(true.strip().split(), test.strip().split())
            
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.test_main']
    unittest.main()
