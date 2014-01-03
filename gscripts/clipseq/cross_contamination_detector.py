'''
Created on Dec 11, 2013

@author: gpratt
'''
import os
from optparse import OptionParser
from itertools import groupby

import pysam


##Need to write up testing code for this

def correlation(bam_1, bam_2, outbam):
    
    """
    
    bam_1: path to indexed bam file
    bam_2: path to indexed bam file
    
    returns number of matched reads between first and second bam file
    and total number of reads in the first bam file
    """


    total_count = 0
    matched_count = 0
    name1 = os.path.basename(".".join(bam_1.split(".")[:2]))
    name2 = os.path.basename(".".join(bam_2.split(".")[:2]))
    
    with pysam.Samfile(bam_1) as bam_1, pysam.Samfile(bam_2) as bam_2:
        outbam = pysam.Samfile(outbam, 'wh', bam_1)
        
        for read in bam_1:
            total_count += 1
            read_start = read.positions[-1] if read.is_reverse else read.positions[0]
            
            fetched_reads = list(bam_2.fetch(bam_1.getrname(read.tid), 
                                             read_start, 
                                             read_start + 1))
            for fetched_read in fetched_reads:
                
                fetched_start = fetched_read.positions[-1] if fetched_read.is_reverse else fetched_read.positions[0]

                if read.qname.split(":")[0] == fetched_read.qname.split(":")[0] and read_start == fetched_start:
                    matched_count += 1
                    outbam.write(read)
                    break
    outbam.close()
    return matched_count, total_count

if __name__ == '__main__':
    usage = """  detects cross contamination between two samples demultiplexed 
    via the demultiplex barcoded fastq script, after alignment """
    parser = OptionParser(usage)
    parser.add_option("-f", "--bam_1", dest="bam_1", 
                      help="first (barcoded) bam file to look for cross contamination with")
    parser.add_option("-b", "--bam_2", dest="bam_2", 
                      help="second (barcoded) bam file to look for cross contamination with")
    parser.add_option("-o", "--out_file", dest="out_file")

    (options, args) = parser.parse_args()
    outbam = pysam.Samfile(os.path.splitext(options.out_file)[0] + ".sam", "wh", pysam.Samfile(options.bam_1))
    matched_count, total_count = correlation(options.bam_1, options.bam_2, outbam)
    outbam.close()
    name1 = os.path.basename(".".join(options.bam_1.split(".")[:2]))
    name2 = os.path.basename(".".join(options.bam_2.split(".")[:2]))
    with open(os.path.join(options.out_file), 'w') as outfile:
        outfile.write("\t".join(map(str, [name1, name2, matched_count, total_count])) + "\n")
        
