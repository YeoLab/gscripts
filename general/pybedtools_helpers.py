"""
List of functions to help help with pybedtools or improve on its features, eventually should be included into pybedtool package

"""


from collections import defaultdict

import numpy as np
import pybedtools

def small_peaks(feature):
    """

    feature - pybedtools feature

    returns center of clipper called peak (the middle of the narrow start / stop)
    
    """
    feature.start = (int(feature[6]) + int(feature[7])) /  2
    feature.stop = (int(feature[6]) + int(feature[7])) /  2
    feature.name = feature.name.split("_")[0]
    return feature

def get_five_prime_end(feature):

    """

    gets 5' end interval in a strand intelgent way (pybedtools implementation of this doesn't work)

    """

    if feature.strand == "+":
        feature.stop = feature.start
    else:
        feature.start = feature.stop
    return feature

def get_three_prime_end(feature):

    """

    gets 3' end interval in a strand intelgent way (pybedtools implementation of this doesn't work)

    """

    if feature.strand == "+":
        feature.start = feature.stop
    else:
        feature.stop = feature.start
    return feature


def adjust_after_shuffle(interval):

    """

    Adjusts name and strand to correct name and strand after shuffling (assumes use of shuffle_transcriptome method)
    
    """
    #Adjusts name and strand in one to name and strand that was in two
    interval.name = interval[11]
    interval.strand = interval[13]
    return interval


def shuffle_and_adjust(bedtool, incl):

    """
    
    bedtool: bedtool to shuffle
    incl: bedtool to include in
    
    Shuffles bedtool and re-adjusts name and strand of interval to match to new location.
    
    """
    
    
    already_exists = {}
    shuffled_tool = bedtool.shuffle(g="/nas3/yeolab/Genome/ucsc/hg19/hg19.chrom.sizes", incl=incl.fn).intersect(incl, wo=True)
    for interval in shuffled_tool:
        existance_tuple = (interval.chrom, interval.start, interval.start, interval.name)
        
        if existance_tuple not in already_exists:
            already_exists[existance_tuple] = interval
            
    return pybedtools.BedTool(already_exists.values()).each(adjust_after_shuffle).saveas()
                                                                
def closest_by_feature(bedtool, closest_feature):

    """
    
    Returns distance from nearest feature assuming distance is centered on the feature
    
    bedtool - a bedtool to find closest feature of
    closest_feature - bedtool of features to find closest things of
    
    Returns closest bed objects 
    
    Assumes both the bedtools object and the feature are 1bp long so we get the distance from both from their start sites
    """
    
    #feature_dict = {feature.name : feature for feature in closest_feature}
    feature_dict = defaultdict(list)
    for feature in closest_feature:
        feature_dict[feature.name].append(feature)
        
    not_included = []
    distances = []
    for interval in bedtool:
        if interval.name not in feature_dict:
            not_included.append(interval.name)
            continue
        
        best_distance = (np.inf, None)
        for feature in feature_dict[interval.name]:
            

            if feature.strand != interval.strand:
                raise ValueError("Strands not identical\nfeature : %sinterval: %s" % (str(feature), str(interval)))
            
            if feature.strand == "+":
                distance = interval.start - feature.start
            else:
                distance = feature.start - interval.start
                
            if abs(distance) < abs(best_distance[0]):
                best_distance = (distance, feature)
        distances.append("\t".join([str(interval).strip(), str(best_distance[1]).strip(), str(best_distance[0])]))

    return pybedtools.BedTool(distances).saveas()
