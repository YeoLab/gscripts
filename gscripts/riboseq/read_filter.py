#!/usr/bin/python

import argparse
import os
import subprocess

import pybedtools


OFFSET = 14

def get_rpf(bam, out_file, length=29):
    """

    Given a bam file and a list of regions returns a dataframe with the distance of each read from the closest region

    bam -- bam file
    regions -- bed file, genomic regions to get distance from
    half_window_width -- int, distance around region to record

    """
    with HTSeq.BAM_Reader(bam) as sorted_bam, HTSeq.BAM_Writer.from_BAM_Reader(out_file, sorted_bam ) as bam_writer:
        for almnt in sorted_bam:
            try:
                if almnt.iv.length == length:
                    bam_writer.write(almnt)
            except AttributeError:
                bam_writer.write(almnt)


parser = argparse.ArgumentParser(description="filters bam files to only have reads of a specific length, important for cleaning up low quality riboseq data")
parser.add_argument("--bam", help="bam file to adjust", required=True)
parser.add_argument("--out", help="output file", required=True)
parser.add_argument("--length", help="length of read to keep", required=False, default=29)
args = parser.parse_args()

get_rpf(args.bam, args.out, args.length)
