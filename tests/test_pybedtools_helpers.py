'''
Created on May 29, 2013

@author: gabrielp

'''


import unittest

import pybedtools

from gscripts.general.pybedtools_helpers import closest_by_feature, convert_to_mRNA_position

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
        
    def test_convert_to_mRNA_position_failaure(self):
        
        """
        
        Tests faialaure state if gene is not found within location_dict
        
        """
        return
        tool = pybedtools.create_interval_from_list("ENSMUSG1    50    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -    50    50".split())
        location_dict = {"ENSMUSG2" : {"strand" : "+", "regions" : [(30, 40),
                                                                    (10,20)
                                                                    ] 
                                       }
                         }
        self.assertRaises(KeyError, convert_to_mRNA_position, tool, location_dict)
        
    def test_convert_to_mRNA_position_strand_equality(self):
        
        """
        
        Tests failaure mode if strands aren't equal
        
        """
        return
        tool = pybedtools.create_interval_from_list("ENSMUSG1    50    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -    60    60".split())
        location_dict = {"ENSMUSG1" : {"strand" : "+", "regions" : [(0,100),
                                                                    ] 
                                       }
                         }
        
        self.assertRaises(ValueError, convert_to_mRNA_position, tool, location_dict)
        
    def test_convert_to_mRNA_position_placement(self):
        
        """
        
        Makes sure that the placement within a region or list of regions is correct
        
        """
        
        return
        interval = pybedtools.create_interval_from_list("ENSMUSG1    50    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    +    60    60".split())
        location_dict = {"ENSMUSG1" : {"strand" : "+", "regions" : [(0,100),
                                                                    ] 
                                       }
                         }
        
        correct_tool = pybedtools.create_interval_from_list("ENSMUSG1    50    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    +    60    60".split())
        self.assertEqual(convert_to_mRNA_position(interval, location_dict), correct_tool)
        
        interval = pybedtools.create_interval_from_list("ENSMUSG1    50    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -".split())
        location_dict = {"ENSMUSG1" : {"strand" : "-", "regions" : [(0,100),
                                                                    ] 
                                       }
                         }
        
        #individual_fraction, total_fraction
        correct_tool = pybedtools.create_interval_from_list("ENSMUSG1    40    50    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -".split())
        x =  convert_to_mRNA_position(interval, location_dict)
        print x
        self.assertEqual(x, correct_tool)
        
    def test_convert_to_mRNA_position_placement_split(self):    
        
        """
        
        Makes sure that lists of regions works for both positive and negative strands
        
        """
        
        return
        tool = pybedtools.create_interval_from_list("ENSMUSG1    125    127    ENSMUSG1_1_83;ENSMUSG1_6_83    0    +    125    125".split())
        location_dict = {"ENSMUSG1" : {"strand" : "+", "regions" : [(0, 50),
                                                                    (100, 150),
                                                                    ] 
                                       }
                         }
        
        correct_tool = pybedtools.create_interval_from_list("ENSMUSG1    75    77    ENSMUSG1_1_83;ENSMUSG1_6_83    0    +    125    125".split())
        self.assertEqual(convert_to_mRNA_position(tool, location_dict), correct_tool )
        
        tool = pybedtools.create_interval_from_list("ENSMUSG1    25    27    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -    25    25".split())
        location_dict = {"ENSMUSG1" : {"strand" : "-", "regions" : [(100, 150),
                                                                    (0, 50),
                                                                    ] 
                                       }
                         }
        
        correct_tool = pybedtools.create_interval_from_list("ENSMUSG1    73    75    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -    25    25".split())
        self.assertEqual(convert_to_mRNA_position(tool, location_dict), correct_tool)
    
    def test_convert_to_mRNA_position_fail(self):
        
        """ Various attempts to break RNA position and make sure that error are caught """
        return    
        tool = pybedtools.create_interval_from_list("ENSMUSG1    51    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -    10    10".split())
        location_dict = {"ENSMUSG1" : {"strand" : "-", "regions" : [(100, 150),
                                                                    (25,50),
                                                                    ] 
                                       }
                         }
        
        self.assertEqual(convert_to_mRNA_position(tool, location_dict).chrom, "none")
        
        tool = pybedtools.create_interval_from_list("ENSMUSG1    51    60    ENSMUSG1_1_83;ENSMUSG1_6_83    0    -    175    175".split())
        
        self.assertEqual(convert_to_mRNA_position(tool, location_dict).chrom, "none")
        
        pybedtools.BedTool("chr1    x     y    ")
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()
