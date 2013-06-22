#!/usr/local/bin/python2.7
# encoding: utf-8
'''
rnaseq.singe_RPKM -- shortdesc

rnaseq.singe_RPKM is a description

It defines classes_and_methods

@author:     user_name
        
@copyright:  2013 organization_name. All rights reserved.
        
@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from collections import defaultdict
def main(counts, outfile): # IGNORE:C0111
    
    total_reads = 0
    gene_lengths = defaultdict(int)
    gene_counts = defaultdict(int)
    
    if outfile != sys.stdout:
        outfile = open(outfile, 'w')
    
    if counts != sys.stdin:
        counts = open(counts)
        
    for line in counts:
        
        chrom, start, stop, region_count, gene_count, strand, gene_id, frea = line.strip().split()
        start, stop, region_count, gene_count = int(start), int(stop), float(region_count), float(gene_count)
        
        total_reads += float(region_count)
        gene_lengths[gene_id] += stop - start
        gene_counts[gene_id] = gene_count
        
    outfile.write("gene    flag    RPKM\n")
    for gene_id in gene_counts.keys():
        RPK = gene_counts[gene_id] / (gene_lengths[gene_id] / 1000.0)
        RPKM = RPK / (total_reads / 1000000)
        outfile.write("\t".join(map(str, [gene_id, 0, RPKM, "\n"])))
        
if __name__ == "__main__":
    parser = ArgumentParser(description="Calculates RPKM for genes")
    parser.add_argument("-i", "--input", dest="input", help="input counts file (from count_tags), default stdin", 
                        default=sys.stdin)
    parser.add_argument("-o", "--output", dest="output", help="output file, default stdout", 
                        default=sys.stdout)
    # Process arguments
    args = parser.parse_args()

    sys.exit(main(args.input, args.output))