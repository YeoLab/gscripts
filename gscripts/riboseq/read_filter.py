#!/usr/bin/python

import argparse
import os
import pysam

def get_match_length(almnt):
    """
    Gets length of perfect matches from cigar string
    """
    if almnt.cigar is None:
        return 0
    
    return sum(length for cigar, length in almnt.cigar if cigar == 0)

def get_rpf(bam, out_file, min_length=28, max_length=30):
    """

    Given a bam file and a list of regions returns a dataframe with the distance of each read from the closest region

    bam -- bam file
    regions -- bed file, genomic regions to get distance from
    half_window_width -- int, distance around region to record

    """

    with pysam.Samfile(bam, 'rb') as sorted_bam:
        with pysam.Samfile(out_file, 'wb', template = sorted_bam ) as bam_writer:
            for almnt in sorted_bam:
                if min_length <= get_match_length(almnt) <= max_length:
                    bam_writer.write(almnt)


parser = argparse.ArgumentParser(description="filters bam files to only have reads of a specific length, important for cleaning up low quality riboseq data")
parser.add_argument("--bam", help="bam file to adjust", required=True)
parser.add_argument("--out", help="output file", required=True)
parser.add_argument("--min_length", help="length of read to keep", required=False, default=28, type=int)
parser.add_argument("--max_length", help="length of read to keep", required=False, default=30, type=int)

args = parser.parse_args()
get_rpf(args.bam, args.out, args.min_length, args.max_length)
