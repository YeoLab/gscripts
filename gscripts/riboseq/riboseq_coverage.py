#!/usr/bin/python

import argparse
import os
import subprocess

import pybedtools


OFFSET = 14

def five_prime(interval):

    #There could be an off by one bug here, need to think about it a bit
    if interval.strand == "+":
        interval.start = interval.start + OFFSET
        interval.end = interval.start + 1
    else:
        interval.end = interval.end - OFFSET 
        interval.start = interval.end - 1

    return interval

parser = argparse.ArgumentParser(description="Adjusts bam files (or bed files) to correct position for riboseq, outputs bed file of corrected positions")
parser.add_argument("--bam", help="bam file to adjust", required=True)
parser.add_argument("--out", help="output file (bed format)", required=True)
args = parser.parse_args()

pybedtools.BedTool(args.bam).bam_to_bed(stream=True).each(five_prime).saveas(args.out + "tmp") 
#bedtools sort runs into memory issues, can't use it. 
p = subprocess.Popen("sort -k 1,1 -k 2,2n %s > %s" % (args.out + "tmp", args.out), shell=True)
p.wait()
os.remove(args.out + "tmp")