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

def shuffleBedTool(name, bedtool):
    
    """
    
    shuffles a bed tool keeping the number of elements in distal and proximal introns the same,
    IMPORTANT: saves result as shuffling is time consuming, results must be retreaved using get_shuffled_lists

    file is saved as name_shuffled.bed
    
    """
    
    #exit if the file has already been made, return its handle
    if os.path.exists(name + "_shuffled.bed"):
        return pybedtools.BedTool(name + "_shuffled.bed")
    
    #else run through... making this general is going to be a bitch, stopping for now.  
    try:
        prox_elements = bedtool.intersect("/nas3/gpratt/projects/structure/prox_event_detail.BED", u=True)
        
        #removes elements that overlap with both prox and distal (assigns them to to prox)
        not_prox_elements = bedtool.intersect("/nas3/gpratt/projects/structure/prox_event_detail.BED", v=True)
     
        #gets distal elements
        dist_elements = not_prox_elements.intersect("/nas3/gpratt/projects/structure/dist_event_detail.BED", u=True)
         
        #shuffles seperatly
        shuffled_prox = prox_elements.shuffle(g="/nas3/gpratt/ens/hg19.genome", incl="/nas3/gpratt/projects/structure/prox_event_detail.BED")
        shuffled_dist = dist_elements.shuffle(g="/nas3/gpratt/ens/hg19.genome", incl="/nas3/gpratt/projects/structure/dist_event_detail.BED")
        
        
        #saves result file
        return_value =  shuffled_prox.cat(shuffled_dist, postmerge=False).saveas(name + "_shuffled.bed")
        return return_value
        
    except Exception as e :
        print "error occured, ignoring", e 
    
    
def shuffle_dict(dict_1, dict_2, num_shuffles=10):
    
    """
    
    Shuffles num_shuffles times
    for two dictionaries that you want to compare against each other, shuffles them. 
    returns two di
    """
    shuffled_dict_1 = {}
    shuffled_dict_2 = []
    for x in range(num_shuffles):
        for dataset_name, dataset_element in dict_1.items():
            if dataset_name not in shuffled_dict_1:
                shuffled_dict_1[dataset_name] = []
                
            shuffled_dict_1[dataset_name].append(shuffleBedTool(dataset_name + str(x), dataset_element))
            
        for dataset_name, dataset_element in dict_2.items():
            if dataset_name not in shuffled_dict_2:
                shuffled_dict_2[dataset_name] = []
                
            shuffled_dict_2[dataset_name].append(shuffleBedTool(dataset_name + str(x), dataset_element))
        
    
    return shuffled_dict_1, shuffled_dict_2

def intersect_shuffled_lists(dict_1, dict_2, shuffled_list_1, shuffled_list_2):
    
    """
    
    Returns a dictionary of lists of all the shuffled lists and every possible intersection 
    
    """
    
    shuffled_lists = {"shuffledUltra_by_clip" : [],
                      "clip_by_shuffledUltra" : [],
                      "ultra_by_shuffledClip" : [],
                      "shuffledClip_by_ultra" : [],
                      "shuffledUltra_by_shuffledClip" : [],
                      "shuffledClip_by_shuffledUltra" : []}
    
    for x in range(10):

            current_shuffledUltra_by_clip, current_clip_by_shuffledUltra = get_all_overlaps(shuffled_ultraDict, dict_1)
            current_ultra_by_shuffledClip, current_shuffledClip_by_ultra  = get_all_overlaps(dict_2, shuffled_CLIPDict)
            current_shuffledUltra_by_shuffledClip, current_shuffledClip_by_shuffledUltra = get_all_overlaps(shuffled_ultraDict, shuffled_CLIPDict)
            
            shuffled_lists["shuffledUltra_by_clip"].append(current_shuffledUltra_by_clip)
            shuffled_lists["clip_by_shuffledUltra"].append(current_clip_by_shuffledUltra)
            shuffled_lists["ultra_by_shuffledClip"].append(current_ultra_by_shuffledClip)
            shuffled_lists["shuffledClip_by_ultra"].append(current_shuffledClip_by_ultra)
            shuffled_lists["shuffledUltra_by_shuffledClip"].append(current_shuffledUltra_by_shuffledClip)
            shuffled_lists["shuffledClip_by_shuffledUltra"].append(current_shuffledClip_by_shuffledUltra)
 
    return shuffled_lists
            
def get_z_scores(shuffled_lists, ultra_by_clip, clip_by_ultra):
 
    shuffledUltra_by_clip_z = z_score(ultra_by_clip, mean(shuffled_lists["shuffledUltra_by_clip"]), std(shuffled_lists["shuffledUltra_by_clip"]))
    #print shuffledUltra_by_clip_z
    #print (rv.pdf(shuffledUltra_by_clip_z) < .05)

    ultra_by_shuffledClip_z = z_score(ultra_by_clip, mean(shuffled_lists["ultra_by_shuffledClip"]), std(shuffled_lists["ultra_by_shuffledClip"]))
    #print ultra_by_shuffledClip_z
    #print (rv.pdf(ultra_by_shuffledClip_z) < .05)

    shuffledUltra_by_shuffledClip_z = z_score(ultra_by_clip, mean(shuffled_lists["shuffledUltra_by_shuffledClip"]), std(shuffled_lists["shuffledUltra_by_shuffledClip"]))
    #print shuffledUltra_by_shuffledClip_z
    #print (rv.pdf(shuffledUltra_by_shuffledClip_z) < .05)

    clip_by_shuffledUltra_z = z_score(clip_by_ultra, mean(shuffled_lists["clip_by_shuffledUltra"]), std(shuffled_lists["clip_by_shuffledUltra"]))
    #print clip_by_shuffledUltra_z
    #print (rv.pdf(clip_by_shuffledUltra_z) < .05)

    shuffledClip_by_ultra_z = z_score(clip_by_ultra, mean(shuffled_lists["shuffledClip_by_ultra"]), std(shuffled_lists["shuffledClip_by_ultra"]))
    #print shuffledClip_by_ultra_z
    #print (rv.pdf(shuffledClip_by_ultra_z) < .05)

    shuffledClip_by_shuffledUltra_z = z_score(clip_by_ultra, mean(shuffled_lists["shuffledClip_by_shuffledUltra"]), std(shuffled_lists["shuffledClip_by_shuffledUltra"]))
    #print shuffledClip_by_shuffledUltra_z
    #print (rv.pdf(shuffledClip_by_shuffledUltra_z) < .05)
    return shuffledUltra_by_clip_z, ultra_by_shuffledClip_z, shuffledUltra_by_shuffledClip_z, clip_by_shuffledUltra_z, shuffledClip_by_ultra_z, shuffledClip_by_shuffledUltra_z

