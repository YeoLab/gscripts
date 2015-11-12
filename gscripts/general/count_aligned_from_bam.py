#!/usr/bin/env python
"""Count number of reads mapped to each chromosome"""

import sys
import argparse
from collections import Counter

import pysam

def count_chroms(handle):
    """
    
    Counts chromosome section of a same file, sorts and prints by most frequent
    
    """

    total_counts = Counter()
    fraction_counts = Counter()
    for read in handle:
        ref = handle.getrname(read.tid)

        total_counts[ref] += 1
        fraction_counts[ref] += 1.0 / read.get_tag("NH")
    for name, count in total_counts.most_common():
        print name, count, fraction_counts[name]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Counts Counts Chromosomes in sam file (used generally for getting counts of rRNA elements after filtering them), writes results to stdout')
    parser.add_argument('--input', help='input bam file')

    args = parser.parse_args()

    handle = pysam.Samfile(args.input, 'rb')
    count_chroms(handle)

    sys.exit(0)
