'''
Created on Dec 11, 2013

@author: gpratt
'''

import tests
import unittest
from gscripts.clipseq.cross_contamination_detector import correlation

class Test(unittest.TestCase):


    def test_randomer_off_by_one_pos(self):
        """
        test off by one bug same randomer but start off by one pos
        """
        
        bam1 = tests.get_file("test_cross_contamination/positive1.bam")
        bam2 = tests.get_file("test_cross_contamination/positive_off_by_one.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 0)
    def test_randomer_off_by_one_neg(self):
        """
        test off by one bug same randomer by off by one neg
        """
        
        bam1 = tests.get_file("test_cross_contamination/negative1.bam")
        bam2 = tests.get_file("test_cross_contamination/negative_off_by_one.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 0)
    def test_randomer_mismatch_pos(self):
        """
        same start different randomer pos
        """
        
        bam1 = tests.get_file("test_cross_contamination/positive1.bam")
        bam2 = tests.get_file("test_cross_contamination/positive_mismatch.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 0)
    def test_randomer_mismatch_neg(self):
        """
        same start different randomer neg
        """
        
        bam1 = tests.get_file("test_cross_contamination/negative1.bam")
        bam2 = tests.get_file("test_cross_contamination/negative_mismatch.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 0)
    def test_randomer_match_pos(self):
        """
        same start same randomer pos
        """
        
        bam1 = tests.get_file("test_cross_contamination/positive1.bam")
        bam2 = tests.get_file("test_cross_contamination/positive_match.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 1)
    def test_randomer_match_neg(self):
        """
        same start same randomer neg
        """
        
        bam1 = tests.get_file("test_cross_contamination/negative1.bam")
        bam2 = tests.get_file("test_cross_contamination/negative_match.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 1)
        
    def test_randomer_match_neg(self):
        """
        same start same randomer neg, with the other offset
        """
        
        bam1 = tests.get_file("test_cross_contamination/negative2.bam")
        bam2 = tests.get_file("test_cross_contamination/negative_match.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 1)
    def test_duplicate_pos(self):
        """
        same start / with target having both a matching and not matching randomer at that location
        """
        
        bam1 = tests.get_file("test_cross_contamination/positive1.bam")
        bam2 = tests.get_file("test_cross_contamination/positive_duplicate.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 1)
    def test_duplicate_neg(self):
        """
        same start / with target having both a matching and not matching randomer at that location
        """
        
        bam1 = tests.get_file("test_cross_contamination/negative1.bam")
        bam2 = tests.get_file("test_cross_contamination/negative_duplicate.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 1)
    def test_pos_vs_neg(self):
        """
        tests same read with same barcode but different strand
        """
        
        bam1 = tests.get_file("test_cross_contamination/positive1.bam")
        bam2 = tests.get_file("test_cross_contamination/negative1.bam")
        matched, total = correlation(bam1, bam2)
        self.assertEqual(total, 1)
        self.assertEqual(matched, 0)
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()