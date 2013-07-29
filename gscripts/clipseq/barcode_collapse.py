<<<<<<< HEAD
""""

barcode_collapse.py  read in a .bam file where the 
first 9 nt of the read name 
are the barcode and merge reads mapped to the same position that have the same barcode

"""

from collections import Counter
import gzip
from optparse import OptionParser
import sys

import numpy as np
import pysam

def barcode_collapse(inBam, outBam, barcoded, use_stop):
    
    """
    
    Removes reads with same start and same barcode
    
    inBam : location of input bam file
    outBam : location of output bam file
    
    returns two dicts, total_count, a dict of all reads and their assocated barcodes
                       removed_count, a dict of all reads that have been removed and their attached barcodes 
                                       {barcode : count} 
    """
    
    outTotal = open(outBam + ".total.wiggle", 'w')
    outBarcodes = open(outBam + ".barcodes.wiggle", 'w')
    outEntropy = open(outBam + ".entropy.wiggle", 'w')
    inBam = pysam.Samfile(inBam, 'rb')
    
    outBam = pysam.Samfile(outBam, 'wb', template=inBam)
    

    cur_count = 0
    prev_chrom = None

    if use_stop:
        prev_pos = (0, 0)
    else:
        prev_pos = 0

    barcode_set = Counter()
    removed_count = Counter()
    total_count = Counter()
    
    for i, read in enumerate(inBam.fetch()):
        cur_chrom = read.rname
        cur_count += 1
        #paramater options to allow for start and stop and barcodes vs normal
        if use_stop:
            cur_pos = (read.positions[0], read.positions[-1])
        else:
            cur_pos = read.positions[0]
        
        if barcoded:
            barcode = read.qname[:9]
        else:
            barcode = "total"

        #if we advance a position, reset barcode counting
        if not cur_pos == prev_pos:
            #make daddy a wiggle track!
            if prev_chrom != cur_chrom:
                var_step = "variableStep chrom=%s\n" % (inBam.getrname(cur_chrom))
                outTotal.write(var_step)
                outBarcodes.write(var_step)

            out_pos = prev_pos[0] if use_stop else prev_pos
            barcode_counts = np.array(barcode_set.values())
            barcode_probablity = barcode_counts / float(sum(barcode_counts))
            entropy = -1 * sum(barcode_probablity * np.log2(barcode_probablity))
            
            outEntropy.write("\t".join(map(str, [out_pos, entropy])) + "\n")
            outTotal.write("\t".join(map(str, [out_pos, cur_count])) + "\n")
            outBarcodes.write("\t".join(map(str, [out_pos, len(barcode_set)])) + "\n")
            barcode_set = Counter()
            cur_count = 0 

        
        total_count[barcode] += 1
        if barcode in barcode_set:
            removed_count[barcode] += 1
        else:
            outBam.write(read)

            
        barcode_set[barcode] += 1   
        prev_pos = cur_pos
        prev_chrom = cur_chrom

    out_pos = prev_pos[0] if use_stop else prev_pos
    outTotal.write("\t".join(map(str, [out_pos, cur_count])) + "\n")
    outBarcodes.write("\t".join(map(str, [out_pos, len(barcode_set)])) + "\n")
    
    inBam.close()
    outBam.close()
    
    outEntropy.close()
    outTotal.close()
    outBarcodes.close()
    return total_count, removed_count

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bam", dest="bam", help="bam file to barcode collapse")
    parser.add_option("-c", "--randomer", action="store_true", dest="barcoded", help="bam files are iclip barcoded")
    parser.add_option("-s", "--use_stop", action="store_true", dest="use_stop", help="use stop as well as start to figure collapse")
    parser.add_option("-o", "--out_file", dest="out_file")
    parser.add_option("-m", "--metrics_file", dest="metrics_file")
    
    
    (options,args) = parser.parse_args()
    
    if not (options.bam.endswith(".bam")):
        raise TypeError("%s, not bam file" % (options.bam))

    pysam.index(options.bam)
    total_count, removed_count = barcode_collapse(options.bam, 
                                                  options.out_file,
                                                  options.barcoded,
                                                  options.use_stop)

    
    with open(options.metrics_file, 'w') as metrics:
        metrics.write("\t".join(["randomer", "total_count", "removed_count"]) + "\n")
        for barcode in total_count.keys():
            metrics.write("\t".join(map(str, [barcode, total_count[barcode], removed_count[barcode]])) + "\n")

    pysam.index(options.out_file)
    sys.exit(0)
