'''
Created on Dec 11, 2013

@author: gpratt
'''
from collections import defaultdict, Counter
from itertools import groupby, permutations
from functools import partial
from optparse import OptionParser
import os

import HTSeq
import pandas as pd
import pysam

import tempfile
import pysam
import pybedtools
from collections import defaultdict, Counter
from itertools import permutations
import pandas as pd

def editor(read):
    if read.strand == "-":
        read.start = read.stop - 1
    else:
        read.stop = read.start + 1
    return read

#assumes properly sorted, barcode collapsed reads, otherwise everything breaks
def tags_to_dict(tags):
    """

    Converts the tag set to a dictionary

    """
    return {key: value for key, value in tags}


def append_read_group(base, read_group):
    try:
        base.add(read_group)
        return base
    except:
        return set([read_group])


def bedtools_count_contamination(bam_file):
    #Hopefully a faster way of doing this...
    handle, fn = tempfile.mkstemp()
    with pysam.AlignmentFile(bam_file, 'rb') as reads:
        with pysam.AlignmentFile(fn, 'wb', template=reads) as writer:
            for x, read in enumerate(reads):
                read.qname = read.qname + ":" + read.tags[-1][1]
                writer.write(read)

    bedtool = pybedtools.BedTool(fn).bam_to_bed()
    bedtool = bedtool.each(editor).sort().saveas()

    #can be a bit sloppy because I'm processing base by base
    #should build in strandedness, but this should be good enough, probablbity is still low...
    total = Counter()
    combinations = defaultdict(Counter)

    prev_base = 0
    base_dict = defaultdict(list)
    for x, interval in enumerate(bedtool):
        #reset counter, and process
        if interval.start != prev_base:
            for randomer, read_groups in base_dict.items():
                for rg1, rg2 in permutations(read_groups, 2):
                    combinations[rg1][rg2] += 1
            #processing code goes here
            base_dict = defaultdict(list)

        name = interval.name.split(":")

        randomer, read_group = name[0], name[-1]
        base_dict[randomer].append(read_group)

        total[read_group] += 1
        prev_base = interval.start

    return pd.DataFrame(combinations), pd.Series(total)


def genome_count_contamination(bam_file):
    #Load all values into dict
    reads_at_location = defaultdict(lambda: HTSeq.GenomicArray("auto", typecode="O", stranded=True))
    reads = HTSeq.BAM_Reader(bam_file)
    for read in reads:
        try:
            randomer = read.read.name.split(":")[0]
            read_group = tags_to_dict(read.optional_fields)['RG']
            reads_at_location[randomer][read.iv].apply(partial(append_read_group, read_group=read_group))
        except AttributeError as e:
            pass

    #Get them out
    combinations = defaultdict(Counter)
    total = Counter()
    for randomer in reads_at_location.keys():
        #print reads_at_location[randomer]
        for genomic_interval, read_groups in reads_at_location[randomer].steps():
            if not read_groups:
                continue

            for read_group in read_groups:
                total[read_group] += 1
            for rg1, rg2 in permutations(read_groups, 2):
                    combinations[rg1][rg2] += 1

    return pd.DataFrame(combinations), pd.Series(total)

def mark_overlap_for_base(reads):

    """

    Given a iterable of reads returns a boolean matrix of read groups by barcodes of what readgroup has what barcode

    """

    counts = defaultdict(dict)
    for read in reads:
        randomer = read.qname.split(":")[0]
        read_group = tags_to_dict(read.tags)["RG"]
        counts[randomer][read_group] = True
    return pd.DataFrame(counts).fillna(False)

def reads_starting_at_location(reads, loc):
    """

    given list of reads returns all reads starting a given loc (sepereated by positive and negative strand that start at a given location)

    """
    pos_reads = []
    neg_reads = []
    for read in reads:
        read_start = read.positions[-1] if read.is_reverse else read.positions[0]
        if read_start == loc:
            if read.is_reverse:
                neg_reads.append(read)
            else:
                pos_reads.append(read)
    return pos_reads, neg_reads

def count_contamination(bam_file):
    """

    Given bam file with readgroups and barcodes counts numboer of overlapping barcodes per readgroup
    return dataframe of overlaps and series of total counts

    """
    combinations = defaultdict(Counter)
    total = Counter()
    with pysam.Samfile(bam_file) as bam_group:
        #pipeup ignores multimapped reads by default, be careful
        for base in bam_group.pileup():
            #given a base, get all reads that start there, either on the positive or negative strands
            read_start = base.pos
            reads = (read.alignment for read in base.pileups)

            pos_reads, neg_reads = reads_starting_at_location(reads, read_start)
            for read in pos_reads:
                read_group = tags_to_dict(read.tags)["RG"]
                total[read_group] += 1

            for read in neg_reads:
                read_group = tags_to_dict(read.tags)["RG"]
                total[read_group] += 1

            pos_overlap = mark_overlap_for_base(pos_reads)
            neg_overlap =  mark_overlap_for_base(neg_reads)

            #Count both the negative and positive overlap
            for randomer in pos_overlap.columns:
                for rg1, rg2 in permutations(pos_overlap[pos_overlap[randomer]].index, 2):
                    combinations[rg1][rg2] += 1

            for randomer in neg_overlap.columns:
                for rg1, rg2 in permutations(neg_overlap[neg_overlap[randomer]].index, 2):
                    combinations[rg1][rg2] += 1
    return pd.DataFrame(combinations), pd.Series(total)

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
