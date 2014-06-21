__author__ = 'gpratt'

import HTSeq
import pandas as pd


def bed_to_genomic_interval(bed):
    """

    Converts bed file to genomic interval (htseq format) file

    """

    for interval in bed:
        yield HTSeq.GenomicPosition(interval.chrom, interval.start, interval.strand)


def get_closest(bam, regions, half_window_width = 50):
    """

    Given a bam file and a list of regions returns a dataframe with the distance of each read from the closest region

    bam -- bam file
    regions -- bed file, genomic regions to get distance from
    half_window_width -- int, distance around region to record

    """
    profile = {}
    sorted_bam = HTSeq.BAM_Reader(bam.fn)
    for x, tss in enumerate(bed_to_genomic_interval(regions)):
        window = HTSeq.GenomicInterval(tss.chrom, tss.pos - half_window_width, tss.pos + half_window_width, tss.strand)
        for almnt in sorted_bam[window]:
            if almnt.iv.strand == "+":
                read_loc = almnt.iv.start - tss.pos
            else:
                read_loc = tss.pos - almnt.iv.end
            profile[almnt.read.name] = {'dist': read_loc, "length": almnt.iv.length}
    return pd.DataFrame(profile)