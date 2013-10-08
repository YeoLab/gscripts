'''
Created on Jun 21, 2013

@author: gabrielp
'''
import unittest

import tests
from gscripts.rnaseq import count_tags

class Test(unittest.TestCase):

    def test_count_to_regions(self):
        
        """
        
        Tests annotation building
        
        """
        
        result = count_tags.count_to_regions(tests.get_file("count_tags_annotation.bed"))
        
        true_result = {"ENSG1" : {'start' : 1, 'stop' : 300, 
                                  'chrom' : 'chr1', 'strand' : "+", 
                                  "frea" : "0", 'gene_id' : 'ENSG1' ,
                                  "regions" : [(1, 100), (200, 300)]},
                       "ENSG2" : {'start' : 400, 'stop' : 700, 
                                  'chrom' : 'chr1', 'strand' : "-", 
                                  "frea" : "0", 'gene_id' : 'ENSG2',
                                  "regions" : [(600,700), (400, 500)] }}

        self.assertDictEqual(true_result, dict(result))

    def test_count_gene(self):
        
        """
        
        Tests count_gene, makes sure it runs and outputs proper counts
        """
        result = count_tags.count_gene(tests.get_file("test.bam"), {'start' : 1, 'stop' : 500, 
                                  'chrom' : 'chr1', 'strand' : "+", 
                                  "frea" : "0", 'raw_count' : 0, 'gene_id' : 'ENSG1' ,
                                  "regions" : [(1, 100), (399, 500)]}, "none")
                
        self.assertEqual(result[0][0], "ENSG1:1-100")
        self.assertAlmostEqual(result[0][1]['counts'].gene_count, 12, delta=3)
        self.assertAlmostEqual(result[0][1]['counts'].region_count, 4, delta=3)
        self.assertEqual(result[0][1]['start'], 1)
        self.assertEqual(result[0][1]['stop'], 100)
        
        self.assertEqual(result[1][0], "ENSG1:399-500")
        self.assertAlmostEqual(result[1][1]['counts'].gene_count, 12, delta=3)
        self.assertAlmostEqual(result[1][1]['counts'].region_count, 8, delta=3)
        self.assertEqual(result[1][1]['start'], 399)
        self.assertEqual(result[1][1]['stop'], 500)
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
