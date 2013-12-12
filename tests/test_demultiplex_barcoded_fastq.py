'''
Created on Jun 21, 2013

@author: gabrielp
'''
import unittest

import tests
from gscripts.clipseq.demultiplex_barcoded_fastq import reformat_read

class Test(unittest.TestCase):
    def test_reformat_read(self):
        
        """
        Tests basic reformatting of read to output known correct result
        """
        
        barcodes = {"GTTG" : "R01"}

        name = "@M01356\n"
        seq = "CCGGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        plus = "+\n"
        quality = "BBBBCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n"
        barcode, randomer, read = reformat_read(name, seq, plus, quality, barcodes)
        
        self.assertEqual("GTTG", barcode)
        self.assertEqual("CCGGTTGTA", randomer)
        self.assertEqual("@CCGGTTGTA:M01356\nTTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n+\nFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n", read)
        
    def test_variable_length_barcodes(self):
        """
        
        Lets see if a length 6 barcode works
        
        """
        
        barcodes = {"GTTTTG" : "R01"}

        name = "@M01356\n"
        seq = "CCGGTTTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        plus = "+\n"
        quality = "FFBBBBCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n"
        barcode, randomer, read = reformat_read(name, seq, plus, quality, barcodes)
        
        self.assertEqual("GTTTTG", barcode)
        self.assertEqual("CCGGTTTTGTA", randomer)
        self.assertEqual("@CCGGTTTTGTA:M01356\nTTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n+\nFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n", read)
    
    def test_multiple_barcodes(self):
        """
        
        Tests if a length 4 and 6 barcode can play together
        
        """
        
        barcodes = {"GTTG"   : "R00",
                    "GTTTTG" : "R01"}

        name = "@M01356\n"
        seq = "CCGGTTTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        plus = "+\n"
        quality = "FFBBBBCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n"
        barcode, randomer, read = reformat_read(name, seq, plus, quality, barcodes)
        
        self.assertEqual("GTTTTG", barcode)
        self.assertEqual("CCGGTTTTGTA", randomer)
        self.assertEqual("@CCGGTTTTGTA:M01356\nTTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n+\nFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n", read)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    