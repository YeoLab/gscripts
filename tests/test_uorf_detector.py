'''
Created on Jan 21, 2014

@author: gpratt
'''
import unittest

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pybedtools


from gscripts.riboseq.uorf_detector import UORF_detector
import tests

class Test(unittest.TestCase):


    def test_get_five_prime_utr_sequences(self):
        """
        
        Tests 5' UTR dictionary creation
        
        """
        
        gff_file = pybedtools.BedTool(tests.get_file("gff_uorf_test.gff"))
        fa_file = tests.get_file("gff_uorf_test.fa")
        intervals = pybedtools.BedTool("""  chr1    5    10    ENSG1    0    +
                                chr1    15    20    ENSG2    0    -
                                chr1    25    30    ENSG3    0    +
                                chr1    35    40    ENSG3    0    +
                                chr1    45    50    ENSG4    0    -
                                chr1    55    60    ENSG4    0    -

                              """, from_string=True)
        detector = UORF_detector()
        test_dict = detector.get_five_prime_utr_sequences(gff_file, fa_file)
        test_dict = {name : [(interval, str(seq.seq)) for interval, seq in tuple] 
                     for name, tuple in test_dict.items()}
        true_dict = {
                       "ENSG1" : [(intervals[0], "GGGGG")],
                       "ENSG2" : [(intervals[1], "AAAAA")],
                       "ENSG3" : [(intervals[2], "GGGGG"), (intervals[3], "TTTTT")],
                       "ENSG4" : [(intervals[5], "GAAAA"), (intervals[4], "CCCCG")],
                       }
        
        self.assertDictEqual(test_dict, true_dict)
    #Test single 5' UTR
    #Test strung together 5' UTR
    #Bonus 5' UTR with exon for finding of 5' UTRs that go beyond start codon
    
    def test_get_uorf_start_stop_baisc(self):
        
        
        """
        
        Test 3 reading frames orf detection
        
        """
        
        intervals = pybedtools.BedTool("""  chr1    1    9    ENSG1    0    +
                                        """, from_string=True)
        
        test_dict = {
                       "ENSG1" : [(intervals[0], SeqRecord(Seq("ATGGGGTAG", IUPAC.unambiguous_dna)))],
                       "ENSG2" : [(intervals[0], SeqRecord(Seq("AATGGGGTAG", IUPAC.unambiguous_dna)))],
                       "ENSG3" : [(intervals[0], SeqRecord(Seq("AAATGGGGTAG", IUPAC.unambiguous_dna)))],
                       }
        
        detector = UORF_detector()
        test_result = detector.get_uorf_start_stop(test_dict, uorf_length=0)
        
        true_result = pybedtools.BedTool(
            [
                 pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "1", "3", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_start:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "7", "9", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_end:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "2", "4", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_start:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "8", "10", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_end:ENSG1;Parent=ENSG1" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "3", "5", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_start:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "9", "11", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_end:ENSG1;Parent=ENSG1" ]),
             ]
                           )
        true_result = [str(x) for x in true_result]
        test_result = [str(x) for x in test_result]
        self.assertListEqual(true_result, test_result)
        
    def test_get_uorf_start_stop_baisc_negative(self):
        
        
        """
        
        Test 3 reading frames orf detection, negative strand
        
        """
        
        intervals = pybedtools.BedTool("""  chr1    1    9    ENSG1    0    -
                                            chr1    1    10    ENSG2    0    -
                                            chr1    1    11    ENSG3    0    -
                                        """, from_string=True)
        
        test_dict = {
                       "ENSG1" : [(intervals[0], SeqRecord(Seq("ATGGGGTAG", IUPAC.unambiguous_dna)))],
                       "ENSG2" : [(intervals[1], SeqRecord(Seq("AATGGGGTAG", IUPAC.unambiguous_dna)))],
                       "ENSG3" : [(intervals[2], SeqRecord(Seq("AAATGGGGTAG", IUPAC.unambiguous_dna)))],
                       }
        detector = UORF_detector()
        test_result = detector.get_uorf_start_stop(test_dict, uorf_length=0)
        
        true_result = pybedtools.BedTool(
            [
                 pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "7", "9", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_start:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "1", "3", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_end:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "7", "9", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_start:ENSG2;Parent=ENSG2" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "1", "3", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_end:ENSG2;Parent=ENSG2" ]),
             
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "7", "9", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_start:ENSG3;Parent=ENSG3" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "1", "3", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_end:ENSG3;Parent=ENSG3" ]),
             ]
                           )
        test_result = [str(x) for x in test_result]
        true_result = [str(x) for x in true_result]
        
        self.assertListEqual(test_result, true_result)

        
    def test_get_uorf_start_stop_split(self):
        
        
        """
        
        start / stop codon detection spanning introns
        
        """
        
        intervals = pybedtools.BedTool("""  chr1    1    6    ENSG1    0    +
                                            chr1    12    18    ENSG1    0    +
                                            chr1    1    6    ENSG2    0    +
                                            chr1    12    19    ENSG2    0    +
                                            chr1    1    6    ENSG3    0    +
                                            chr1    12    20    ENSG3    0    +
                                        """, from_string=True)
        
        test_dict = {
                       "ENSG1" : [(intervals[0], SeqRecord(Seq("GGGATG", IUPAC.unambiguous_dna))), 
                                  (intervals[1], SeqRecord(Seq("GGGTAG", IUPAC.unambiguous_dna)))],
                       "ENSG2" : [(intervals[2], SeqRecord(Seq("GGATGG", IUPAC.unambiguous_dna))),
                                  (intervals[3], SeqRecord(Seq("GGTAGG", IUPAC.unambiguous_dna)))],
                       "ENSG3" : [(intervals[4], SeqRecord(Seq("GATGGG", IUPAC.unambiguous_dna))),
                                  (intervals[5], SeqRecord(Seq("GTAGGG", IUPAC.unambiguous_dna)))],
                       }
        
        detector = UORF_detector()
        test_result = detector.get_uorf_start_stop(test_dict, uorf_length=0)
        
        true_result = pybedtools.BedTool(
            [
                 pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "4", "6", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_start:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "15", "17", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_end:ENSG1;Parent=ENSG1" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "3", "5", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_start:ENSG2;Parent=ENSG2" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "14", "16", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_end:ENSG2;Parent=ENSG2" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "2", "4", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_start:ENSG3;Parent=ENSG3" ]),
                
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "13", "15", ".", 
                                                      "+", ".", 
                                                      "ID=uORF_end:ENSG3;Parent=ENSG3" ]),
             ]
                           )
        true_result = [str(interval) for interval in true_result]
        test_result = [str(interval) for interval in test_result]
        print test_result
        print true_result
        self.assertEqual(true_result, test_result)

    def test_get_uorf_start_stop_split_negative(self):
        
        
        """
        
        start / stop codon detection spanning introns negative strand
        
        """
        
        intervals = pybedtools.BedTool("""  chr1    12    18    ENSG1    0    -
                                            chr1    1    6    ENSG1    0    -
                                            chr1    12    18    ENSG2    0    -
                                            chr1    1    6    ENSG2    0    -
                                            chr1    12    18    ENSG3    0    -
                                            chr1    1    6    ENSG3    0    -
                                        """, from_string=True)
        
        test_dict = {
                       "ENSG1" : [(intervals[0], SeqRecord(Seq("GGGATG", IUPAC.unambiguous_dna))), 
                                  (intervals[1], SeqRecord(Seq("GGGTAG", IUPAC.unambiguous_dna)))],
                       "ENSG2" : [(intervals[2], SeqRecord(Seq("GATGGG", IUPAC.unambiguous_dna))),
                                  (intervals[3], SeqRecord(Seq("GTAGGG", IUPAC.unambiguous_dna)))],
                       "ENSG3" : [(intervals[4], SeqRecord(Seq("GGATGG", IUPAC.unambiguous_dna))),
                                  (intervals[5], SeqRecord(Seq("GGTAGG", IUPAC.unambiguous_dna)))]
                       }
        
        detector = UORF_detector()
        test_result = detector.get_uorf_start_stop(test_dict, uorf_length=0)
        
        true_result = pybedtools.BedTool(
            [
                 pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "13", "15", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_start:ENSG1;Parent=ENSG1" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "1", "3", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_end:ENSG1;Parent=ENSG1" ]),
             
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "15", "17", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_start:ENSG2;Parent=ENSG2" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "3", "5", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_end:ENSG2;Parent=ENSG2" ]),
             
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_start", "14", "16", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_start:ENSG3;Parent=ENSG3" ]),
                pybedtools.create_interval_from_list(["chr1", "protein_coding", 
                                                      "uORF_end", "2", "4", ".", 
                                                      "-", ".", 
                                                      "ID=uORF_end:ENSG3;Parent=ENSG3" ]),
             ]
                           )
        true_result = [str(interval) for interval in true_result]
        test_result = [str(interval) for interval in test_result]
        print test_result
        print true_result
        self.assertEqual(test_result, true_result)
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
    
    
    for five_prime_utr, sequence in zip(filtered_UTR5, sequences):
        five_prime_utr_dict[five_prime_utr.name].append((five_prime_utr, sequence))
    
