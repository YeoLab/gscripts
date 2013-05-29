'''
Created on May 29, 2013

@author: gabrielp
'''


import unittest

import pybedtools

from general.pybedtools_helpers import closest_by_feature

class Test(unittest.TestCase):


    def test_closest_by_feature_positive_strand(self):
        
        """
        
        Tests custom closest feature function, need 4 tests, + on tool - on
        features and all other combinations
        
        """
        
        #makes bedtool to test +
        bedtool = pybedtools.BedTool("chr1 100 200 foo 0 +", from_string=True)
        
        closest_feature = pybedtools.BedTool("""chr1 111 200 foo 0 +\n
                                                chr1 90  150 foo 0 +""", from_string=True)
        
        result = closest_by_feature(bedtool, closest_feature)
        
        self.assertEqual("\t".join("chr1 100 200 foo 0 + chr1 90  150 foo 0 + 10".split()) + "\n", str(result))


    def test_closest_by_feature_negative_strand(self):
        
        """
        
        Tests negative strand on closets_by_feature 
        
        """
                    #makes bedtool to test -
        bedtool = pybedtools.BedTool("chr1 100 200 foo 0 -", from_string=True)
        
        closest_feature = pybedtools.BedTool("""chr1 111 200 foo 0 -\n
                                                chr1 90  150 foo 0 -""", from_string=True)
        
        result = closest_by_feature(bedtool, closest_feature)
        
        self.assertEqual("\t".join("chr1 100 200 foo 0 - chr1 90  150 foo 0 - -10".split()) + "\n", str(result))

    def test_closest_by_feature_multiple_beds(self):
        
        """
        
        Tests multiple features and multiple beds
        
        """
        
                            #makes bedtool to test -
        bedtool = pybedtools.BedTool("""chr1 299 300 bar 0 -
                                        chr1 301 400 bar 0 -
                                        chr1 89 200 foo 0 +
                                        chr1 91 200 foo 0 +""", from_string=True)
        
        closest_feature = pybedtools.BedTool("""chr1 310 400 bar 0 -\n
                                                chr1 300 350 bar 0 -\n
                                                chr1 111 200 foo 0 +\n
                                                chr1 90  150 foo 0 +\n""", from_string=True)
        
        result = closest_by_feature(bedtool, closest_feature)
        
        self.assertEqual(str(pybedtools.BedTool("""chr1 299 300 bar 0 - chr1 300 350 bar 0 - 1\n
                                                   chr1 301 400 bar 0 - chr1 300 350 bar 0 - -1\n
                                                   chr1 89 200 foo 0 + chr1 90  150 foo 0 + -1\n
                                                   chr1 91 200 foo 0 + chr1 90  150 foo 0 + 1\n""", from_string=True)), str(result))

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()