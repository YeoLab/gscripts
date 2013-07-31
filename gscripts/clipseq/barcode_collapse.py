""""

barcode_collapse.py  read in a .bam file where the 
first 9 nt of the read name 
are the barcode and merge reads mapped to the same position that have the same barcode
<<<<<<< HEAD
=======

"""

from collections import Counter
import gzip
from optparse import OptionParser
import sys

import numpy as np
import pysam

def barcode_collapse(inBam, outBam, barcoded):
    
    """
    
    Removes reads with same start and same barcode
    
    inBam : location of input bam file
    outBam : location of output bam file
    
    returns two dicts, total_count, a dict of all reads and their assocated barcodes
                       removed_count, a dict of all reads that have been removed and their attached barcodes 
                                       {barcode : count} 
    """
    inBam = pysam.Samfile(inBam, 'rb')
    
    outBam = pysam.Samfile(outBam, 'wb', template=inBam)
    
    prev_pos = (0, 0)
    
    #want total number of reads assocated with each barcode
    #number of reads removed assocated with each barcode
    #some sort of histogram measure of where duplicates are being removed from the most
    #ideally this would be a wiggle track, but for now I'll just keep track of a histogram / before and after trick
    barcode_set = set([])
    removed_count = Counter()
    total_count = Counter()
    
    for i, read in enumerate(inBam.fetch()):
        cur_pos = (read.positions[0], read.positions[-1])
        
        if barcoded:
            barcode = read.qname[:9]
        else:
            barcode = "total"
            
        #if we advance a position, reset barcode counting
        if not cur_pos == prev_pos:
            barcode_set = set([])
            

        
        total_count[barcode] += 1
        if barcode in barcode_set:
            removed_count[barcode] += 1
        else:
            outBam.write(read)

            
        barcode_set.add(barcode)    
        prev_pos = cur_pos
    
    inBam.close()
    outBam.close()
    return total_count, removed_count

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bam", dest="bam", help="bam file to barcode collapse")
    parser.add_option("-c", "--barcoded", action="store_true", dest="barcoded", help="bam files are iclip barcoded")
    parser.add_option("-o", "--out_file", dest="out_file")
    parser.add_option("-m", "--metrics_file", dest="metrics_file")
    
    
    (options,args) = parser.parse_args()
    
    if not (options.bam.endswith(".bam")):
        raise TypeError("%s, not bam file" % (options.bam))
    
    total_count, removed_count = barcode_collapse(options.bam, options.out_file, options.barcoded)
    
    with open(options.metrics_file, 'w') as metrics:
        metrics.write("\t".join(["barcode", "total_count", "removed_count"]) + "\n")
        for barcode in total_count.keys():
            metrics.write("\t".join(map(str, [barcode, total_count[barcode], removed_count[barcode]])) + "\n")

    sys.exit(0)

