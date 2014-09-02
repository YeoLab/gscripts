'''
Created on Jul 30, 2013

@author: gabrielp
'''
import unittest

import tests
from gscripts.clipseq import barcode_collapse
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
        total_count, removed_count = barcode_collapse.barcode_collapse(inBam, outBam, True)
        
        
        true_total_count = {"AAGGGTCGC": 3,
                            "AAGGGTCGT": 4}
        
        true_removed_count = {"AAGGGTCGC": 1,
                            "AAGGGTCGT": 2}

        self.assertDictEqual(true_total_count, total_count)
        self.assertDictEqual(true_removed_count, removed_count)
        
    def test_barcode_collapse_not_barcoded(self):
        """
        Tests duplicate removal for non barcoded samples
        """
        
        inBam = tests.get_file("test_barcode_collapse.bam")
        outBam = tests.get_file("test_barcode_collapse.bc.bam")
        total_count, removed_count = barcode_collapse.barcode_collapse(inBam, outBam, False)
        
        
        true_total_count = {"total": 7}
        
        true_removed_count = {"total": 5}

        self.assertDictEqual(true_total_count, total_count)
        self.assertDictEqual(true_removed_count, removed_count)

    def test_hamming(self):
        """
        Tests hamming distance
        :return:
        """
        self.assertEqual(0, barcode_collapse.hamming("AAAC", "AAAC"))
        self.assertEqual(0, barcode_collapse.hamming("AAAN", "AAAC"))

        self.assertEqual(1, barcode_collapse.hamming("AAAA", "AAAC"))

    def test_calculate_p_read_given_barcode(self):
            self.assertEqual(barcode_collapse.calculate_p_read_given_barcode("AAAC", "AAAC", .05), 1)
            self.assertEqual(barcode_collapse.calculate_p_read_given_barcode("AAAC", "AAAG", .05), .05)

    def test_calculate_p_barcode_given_read(self):
        p_read_given_barcode = {"AAAG": {"AAAG": 1.0, "AAAC": .1},
                                "AAAC": {"AAAG": .1, "AAAC": 1.0, },
                                }
        barcodes_frequency = {"AAAG": .1, "AAAC": .9}

        self.assertAlmostEqual(barcode_collapse.calculate_p_barcode_given_read("AAAC", "AAAG", p_read_given_barcode, barcodes_frequency), .47,
                               delta=2)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()