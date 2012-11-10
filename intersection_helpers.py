'''
Created on Nov 9, 2012

@author: gpratt
'''
import pybedtools
import numpy as np
from scipy import stats
from scipy.stats import norm
import random
from pandas import DataFrame, Series
import os
pybedtools.set_tempdir("/nas/nas0/scratch/gpratt/pybedtools_tmp")
import math
from scipy.stats import norm
rv = norm()

def get_all_overlaps(set1, set2):
    
    """
    
    given dicts set1 {name : pybedtool, name : pybedtool...}
                set2 {name : pybedtool, name : pybedtool...}
    
    returns overlapping counts of all pairwise combinations (in two data frames)
    set1_by_set2, set2_by_set1
     
    """
    set1_by_set2 = {}
    set2_by_set1 = {}
    for set1Name, set1Element in set1.items():
        
        set1_lst  = []
        set2_lst = []
        for set2Name, set2Element in set2.items():
            set1_lst.append(len(set1Element.intersect(set2Element, u=True, stream=True)))
            set2_lst.append(len(set2Element.intersect(set1Element, u=True, stream=True)))
        set1_by_set2[set1Name] = Series(set1_lst, index = set2.keys())
        set2_by_set1[set1Name] = Series(set2_lst, index = set2.keys())
    
    return DataFrame(set1_by_set2), DataFrame(set2_by_set1)


def z_score(x, u, s):
    """
    
    calculate z score x - value, u - sample mean, s - std devation
    
    """
    
    return (x - u) / s



def std(array):
    """
    
    vector wise standard deviation calculation 
    sqrt((sum_x2 / n) - (mean * mean)) 
    
    """
    result = []
    for arr in array:
        result.append(pow(arr - np.mean(array),2))
    result = np.sqrt(sum(result) / 10)

    return result

def get_raw_percents(ultraDict, CLIPDict):
    """
    
    given two dicts prints the fraction of one dict that intersects with the other
    
    """
    ultra_by_clip, clip_by_ultra = get_all_overlaps(ultraDict, CLIPDict)
    print "# ultra conserved elements that intersect with CLIP-seq data"
    print ultra_by_clip #devide by ultra conserved counts

    print "# CLIP-seq peaks that intersect with ultra conserved elements"
    print clip_by_ultra #devide by clip counts
    clip_lens = Series([ len(x) * 1.0 for x in CLIPDict.values() ], index = CLIPDict.keys())
    ultra_lens = Series([ len(x) * 1.0 for x in ultraDict.values() ], index= ultraDict.keys()) 

    print "# CLIP-seq peaks"
    print clip_lens

    print "# ultra conserved elements"
    print ultra_lens

    print "% ultra conserved elements that intersect with CLIP-seq data"
    ultra_by_clip_fraction = (ultra_by_clip / ultra_lens) * 100
    print ultra_by_clip_fraction
    print "% CLIP-seq peaks that intersect with ultra conserved elements"
    clip_by_ultra_fraction =  (clip_by_ultra.T / clip_lens).T * 100
    print clip_by_ultra_fraction
    return ultra_by_clip, ultra_by_clip_fraction, clip_by_ultra, clip_by_ultra_fraction
