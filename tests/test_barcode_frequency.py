'''
Created on Dec 11, 2013

@author: gpratt
'''

from collections import defaultdict
import unittest

import pandas as pd
from gscripts.clipseq.barcode_frequency import handle_seq

class Test(unittest.TestCase):
    """
    Tests barcode frequency methods
    """

    def test_handle_seq(self):
        """
        Test tests one read with two barcodes, a 6mer and a 4mer with the same first 4 bases, everything else should
        report 0
        """
        result_dict = defaultdict(dict)
        barcode_df = pd.DataFrame({"ATGCGC" : {"id" : "R01" }, "ATGC" : {"id" : "R02"}}).T
        seq = "ATGCGCAAA"
        handle_seq(seq, barcode_df, result_dict)
        correct = defaultdict(dict)
        correct["ATGCGC"] = {0 : 1} 
        correct["ATGC"]   =  {0 : 1} 

        self.assertDictEqual(correct, result_dict)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()