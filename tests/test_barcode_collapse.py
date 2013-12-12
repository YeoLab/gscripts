'''
Created on Jul 30, 2013

@author: gabrielp
'''
import unittest

import tests
from gscripts.clipseq.barcode_collapse import barcode_collapse

class Test(unittest.TestCase):
    """
    Tests for barcode collaposer
    """

    def test_barcode_collapse_bacoded(self):
        """
        Tests on duplciate removal for barcoded samples
        """
        
        inBam = tests.get_file("test_barcode_collapse.bam")
        outBam = tests.get_file("test_barcode_collapse.bc.bam")
        total_count, removed_count = barcode_collapse(inBam, outBam, True)
        
        
        true_total_count = {"AAGGGTCGC" : 3, 
                            "AAGGGTCGT" : 4}
        
        true_removed_count = {"AAGGGTCGC" : 1, 
                            "AAGGGTCGT" : 2}

        self.assertDictEqual(true_total_count, total_count)
        self.assertDictEqual(true_removed_count, removed_count)
        
    def test_barcode_collapse_not_barcoded(self):
        """
        Tests duplicate removal for non barcoded samples
        """
        
        inBam = tests.get_file("test_barcode_collapse.bam")
        outBam = tests.get_file("test_barcode_collapse.bc.bam")
        total_count, removed_count = barcode_collapse(inBam, outBam, False)
        
        
        true_total_count = {"total" : 7}
        
        true_removed_count = {"total" : 5}

        self.assertDictEqual(true_total_count, total_count)
        self.assertDictEqual(true_removed_count, removed_count)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()