__author__ = 'gpratt'
'''
Created on Jun 21, 2013

@author: gabrielp
'''
import unittest

import tests
from gscripts.clipseq.demux_paired_end import reformat_read, read_has_barcode

class Test(unittest.TestCase):

    def test_read_has_barcode(self):
        """
        Test hamming distance with hamming of 0
        :return:
        """
        barcodes = ['GTTG', 'AAAA']
        seq_1 = "GTTGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        barcode = read_has_barcode(barcodes, seq_1)
        self.assertEqual(barcode, barcodes[0])

    def test_read_has_barcode_hamming_1(self):
        """
        Test hamming distance with hamming of 0 and 1
        :return:
        """

        barcodes = ['GTTG', 'AAAA']
        seq_1 = "ATTGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        barcode = read_has_barcode(barcodes, seq_1, max_hamming_distance=0)
        self.assertEqual(barcode, None)

        seq_1 = "ATTGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        barcode = read_has_barcode(barcodes, seq_1, max_hamming_distance=1)
        self.assertEqual(barcode, barcodes[0])

    def test_variable_length_barcodes(self):
        """
        Test hamming distance with hamming of 0 and 1
        :return:
        """

        barcodes = ['GTTG', 'AAAAA']
        seq_1 = "ATTGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        barcode = read_has_barcode(barcodes, seq_1, max_hamming_distance=1)
        self.assertEqual(barcode, barcodes[0])

        seq_1 = "AAAAATGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        barcode = read_has_barcode(barcodes, seq_1, max_hamming_distance=1)
        self.assertEqual(barcode, barcodes[1])

    def test_reformat_read(self):

        """
        Tests basic reformatting of read to output known correct result
        """

        barcodes = {"GTTG": "R01"}

        name_1 = "@M01356\n"
        seq_1 = "GTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        plus_1 = "+\n"
        quality_1 = "BBBBCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n"

        name_2 = "@M01356\n"
        seq_2 = "CCGGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n"
        plus_2 = "+\n"
        quality_2 = "BBBBCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n"
        barcode, randomer, read_1, read_2 = reformat_read(name_1, seq_1, plus_1, quality_1,
                                                name_2, seq_2, plus_2, quality_2,
                                                barcodes)

        self.assertEqual("GTTG", barcode)
        self.assertEqual("CC", randomer)
        self.assertEqual(read_1,
                         "@CC:M01356\nTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n+\nCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n")
        self.assertEqual(read_2,
                         "@CC:M01356\nGGTTGTATTTCATTCTGCCCAGAGCAAAATACATGTGACAAAA\n+\nBBCBCFFFFFGGGGGGGGGGHHHHHHHHHHHHHHHHHHHHHHH\n")

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
