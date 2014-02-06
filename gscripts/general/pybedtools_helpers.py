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
    feature.start = (int(feature[6]) + int(feature[7])) / 2
    feature.stop = (int(feature[6]) + int(feature[7])) / 2
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
    shuffled_tool = bedtool.shuffle(
        g="/nas3/yeolab/Genome/ucsc/hg19/hg19.chrom.sizes",
        incl=incl.fn).intersect(incl, wo=True)
    for interval in shuffled_tool:
        existance_tuple = (
        interval.chrom, interval.start, interval.start, interval.name)

        if existance_tuple not in already_exists:
            already_exists[existance_tuple] = interval

    return pybedtools.BedTool(already_exists.values()).each(
        adjust_after_shuffle).saveas()


def make_ucsc_chr(interval):
    """

    Converts interval from ENSEMBL chroms to UCSC chroms

    (appends str to each chrom)

    """

    interval.chrom = "chr" + interval.chrom
    return interval


def rename_to_gene_from_dict(interval, transcript_gene_dict):
    """

    Renames interval name field (and second column because are working on gff files) to the parent gene
    Generally takes a dict {transcript : gene} and converts transcript names to gene names
    
    """

    interval[2] = transcript_gene_dict[interval.attrs['Parent']]
    interval.name = transcript_gene_dict[interval.attrs['Parent']]
    return interval


def get_single_gene_name(interval):
    """

    Converts the name file from an interval to just have one name (if multiple names are combined via a ;

    """

    interval.name = interval.name.split(";")[0]
    return interval


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

            #should throw in stronger error checking here, this is due to slightly different gene annotation approaches being used.  
            if feature.strand != interval.strand or feature.chrom != interval.chrom:
                #continue
                print interval.strand, feature.strand
                raise ValueError(
                    "Strands not identical\nfeature : %sinterval: %s" % (
                    str(feature), str(interval)))

            if feature.strand == "+":
                distance = interval.start - feature.start
            else:
                distance = feature.start - interval.start

            if abs(distance) < abs(best_distance[0]):
                best_distance = (distance, feature)
            #avoids problem of skipping sections
        if best_distance[1] is not None:
            distances.append("\t".join(
                [str(interval).strip(), str(best_distance[1]).strip(),
                 str(best_distance[0])]))

    return pybedtools.BedTool(distances).saveas()


def convert_to_mRNA_position(interval, gene_model):
    """
    
    Returns distance from nearest feature assuming distance is centered on the feature
    
    bedtool - a bedtool to find closest feature of
    gene_model - dict of lists of pybedtools 
    
    Returns bed objects mapped to mRNA position instead of genomic position
    
    Assumes both the bedtools object and the feature are 1bp long so we get the distance from both from their start sites
    
    Negative strand gets modified to be positive strand like, this will fuck with closest bed
    need to do everything on the positive strand from here on out
    """

    #feature_dict = {feature.name : feature for feature in closest_feature}


    if interval.chrom not in gene_model:
        interval.chrom = "none"
        return interval
        #raise KeyError(interval.chrom + " not in current as stucture dict ignoring cluster ")

    #gene model - dict of list of intervals, always at least length 1 
    if not interval.strand == gene_model[interval.chrom][0].strand:
        interval.chrom = "none"
        return interval
        #raise ValueError("strands not the same, there is some issue with gene annotations")

    running_length = 0
    exons = pybedtools.BedTool(gene_model[interval.chrom])
    merged_exons = exons.merge(nms=True, scores='max', s=True)
    _ = pybedtools.helpers.close_or_delete(exons)
    sortmerge = merged_exons.sort().saveas()
    _ = pybedtools.helpers.close_or_delete(merged_exons)
    exons = [i for i in sortmerge]
    _ = pybedtools.helpers.close_or_delete(sortmerge)
    if exons[0].strand == '-':
        exons = exons[::-1]
    for region in exons:

        if interval.start >= int(region.start) and interval.start <= int(
                region.stop):
            if interval.strand == "+":
                tmp_start = running_length + (interval.start - region.start)
                tmp_end = running_length + (interval.end - region.start)

            elif interval.strand == "-": #need the temps for swaping start and end
                tmp_start = running_length + (region.stop - interval.end)
                tmp_end = running_length + (region.stop - interval.start)

            if int(tmp_start) == 0 or int(tmp_end) == 0:
                tmp_start = 1
                tmp_end = 1

            interval.start = tmp_start
            interval.end = tmp_end

            return interval
        running_length += region.length
    interval.chrom = "none"
    return interval

def to_bed(x):
    return x.chrom, x.start, x.stop, x.attributes['gene_id'], "0", x.strand
