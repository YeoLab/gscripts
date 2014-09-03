"""

barcode_collapse.py  read in a .bam file where the 
first 9 nt of the read name 
are the barcode and merge reads mapped to the same position that have the same barcode

"""


from collections import Counter, OrderedDict, defaultdict
import itertools
from optparse import OptionParser
import sys

import numpy as np
import pysam


#if we decide to add back in the start and stop functionality again keep in mind that
#reads are sorted by start site only, not start and stop site, so will need to use a pipeup based stragety
#for removing reads with same start and stop, simply iterating will not work

def hamming(word1, word2):
    """
    Gets hamming distance between two words, this is an odd implementation because
    we actually want Ns to appear similar to other barcodes as this makes the results more stringent, not less
    :param word1:
    :param word2:
    :return:
    """
    return sum(a != b and not (a == "N" or b == "N") for a, b in zip(word1, word2))


def calculate_p_read_given_barcode(read, barcode, error_rate):
    """
    Give a read, barcode and error rate predicts the probablity that the given read came from the given barcode
    :param read: string
    :param barcode: string
    :param error_rate: float
    :return:
    """
    return np.power(error_rate, hamming(read, barcode))


def memoize_p_read_given_barcode(barcodes, error_rate):
    result = defaultdict(dict)
    for read, barcode in itertools.product(barcodes, repeat=2):
        result[read][barcode] = calculate_p_read_given_barcode(read, barcode, error_rate)
    return result


def memoize_p_barcode_given_read(barcodes, p_read_given_barcode, barcodes_frequency):
    result = defaultdict(dict)
    for barcode, read in itertools.product(barcodes, repeat=2):
        result[barcode][read] = calculate_p_barcode_given_read(barcode, read, p_read_given_barcode, barcodes_frequency)
    return result


def memoize_barcodes_frequency(barcodes, reads, p_barcode_given_read):
    return {barcode: calculate_barcode_frequency(barcode, reads, p_barcode_given_read) for barcode in barcodes}


def calculate_p_read_given_barcode(read, barcode, error_rate):
    """
    Give a read, barcode and error rate predicts the probablity that the given read came from the given barcode
    :param read: string
    :param barcode: string
    :param error_rate: float
    :return:
    """
    return np.power(error_rate, hamming(read, barcode))


def calculate_p_barcode_given_read(barcode, read, p_read_given_barcode, barcodes_frequency):
    """
    Given a barcode, read, dict of P(read|barcode) and a frequency of barcodes, predicts P(barcode|read)
    :param barcode: string
    :param read: string
    :param p_read_given_barcode: defaultdict(read: {barcode: Probablity}
    :param barcodes_frequency: dict Estimated frequency of read in sample
    :return:
    """
    numerator = (p_read_given_barcode[read][barcode] * barcodes_frequency[barcode])
    demoninator = sum((p_read_given_barcode[read][cur_barcode] * barcodes_frequency[cur_barcode]) for cur_barcode in barcodes_frequency)
    return numerator / demoninator


def calculate_barcode_frequency(barcode, reads, p_barcode_given_read):
    result = sum(p_barcode_given_read[barcode][read] for read in reads)
    return result / len(reads)


def em_collapse_base(reads, outBam, randomer):
    print "collapsing base"
    barcode_set = {}
    barcodes = []
    for read in reads:
        barcode = read.qname.split(":")[0] if randomer else "total"
        barcodes.append(barcode)
        if barcode not in barcode_set:
            barcode_set[barcode] = read

    barcodes_count = Counter(barcodes)

    barcodes_frequency = {barcode: float(count) / sum(barcodes_count.values()) for barcode, count in barcodes_count.items()}
    error_rate = .05
    p_read_given_barcode = memoize_p_read_given_barcode(barcodes_count, error_rate)
    p_barcode_given_read = memoize_p_barcode_given_read(barcodes_count, p_read_given_barcode, barcodes_frequency)

    #check if each tag exists:
    for barcode, bam_read in barcode_set.items():
        q = sum(np.log10(1 - p_barcode_given_read[barcode][read]) for read in reads)

        if q >= 50:
            outBam.write(read)

    return barcodes_count

def collapse_base(reads, outBam, randomer, total_count, removed_count, max_hamming_distance, em=False):
    
    """
    
    Given a list of reads a base collapses them based on barcodes or just collapses them
    NOT STRAND SPECIFIC, assumes all reads in list are from same positions, and collapses purely based on barcode
    
    Oddly can't return and add two counters together, it creates another counter and causes the entire thing to be very slow
    """
    if em:
        em_collapse_base(reads, outBam, randomer)
    else:
        barcode_set = Counter()
        for read in reads:
            barcode = read.qname.split(":")[0] if randomer else "total"

            if barcode in barcode_set:
                removed_count[barcode] += 1
            else:
                add_read = True
                for cur_barcode in barcode_set:
                    if hamming(barcode, cur_barcode) <= max_hamming_distance:
                        add_read = False
                if add_read:
                    outBam.write(read)

            total_count[barcode] += 1
            barcode_set[barcode] += 1
        return barcode_set


def barcode_collapse(inBam, outBam, randomer, max_hamming_distance=2, em=False):
    
    """
    
    Removes reads with same start and same barcode
    
    inBam : location of input bam file
    outBam : location of output bam file
    
    returns two dicts, total_count, a dict of all reads and their assocated barcodes
                       removed_count, a dict of all reads that have been removed and their attached barcodes 
                                       {barcode : count} 
    """
    
    out_total = open(outBam + ".total.wiggle", 'w')
    out_barcodes = open(outBam + ".barcodes.wiggle", 'w')
    out_entropy = open(outBam + ".entropy.wiggle", 'w')

    inBam = pysam.Samfile(inBam, 'rb')
    
    outBam = pysam.Samfile(outBam, 'wb', template=inBam)
    

    cur_count = 0
    prev_chrom = None
    prev_pos = 0

    #dictionary for handeling reads coming from negative strand
    neg_dict = OrderedDict()
    pos_list = [] #positive we iterate one base at a time, so no need for a dict, a list will do.  
    barcode_set = ([])
    removed_count = Counter()
    total_count = Counter()

    for i, read in enumerate(inBam.fetch()):
        cur_chrom = read.rname
        cur_count += 1
    
        #paramater options to allow for start and stop and barcodes vs normal
        start = read.positions[-1] if read.is_reverse else read.positions[0]
        cur_pos = read.positions[0]

        #if we advance a position, reset barcode counting
        if not cur_pos == prev_pos:
            
            for (key_chrom, key_pos), reads in neg_dict.items():
                if key_pos >= cur_pos:
                    break

                collapse_base(reads, outBam, randomer, total_count, removed_count, max_hamming_distance, em)
                del neg_dict[(key_chrom, key_pos)]


            if prev_chrom != cur_chrom:
                var_step = "variableStep chrom=%s\n" % (inBam.getrname(cur_chrom))
                out_total.write(var_step)
                out_barcodes.write(var_step)
                
                for x in neg_dict.keys():
                    
                    collapse_base(neg_dict[x], outBam, randomer, total_count, removed_count, max_hamming_distance, em)
                    del neg_dict[x]

                assert len(neg_dict) == 0
                
            barcode_set = collapse_base(pos_list, outBam, randomer, total_count, removed_count, max_hamming_distance, em)

            barcode_counts = np.array(barcode_set.values())
            barcode_probablity = barcode_counts / float(sum(barcode_counts))
            entropy = -1 * sum(barcode_probablity * np.log2(barcode_probablity))
            
            out_pos = prev_pos
            out_entropy.write("\t".join(map(str, [out_pos, entropy])) + "\n")
            out_total.write("\t".join(map(str, [out_pos, cur_count])) + "\n")
            out_barcodes.write("\t".join(map(str, [out_pos, len(barcode_set)])) + "\n")
            
            pos_list = []
            cur_count = 0 
        
        if read.is_reverse:
            try:
                neg_dict[(cur_chrom, start)].append(read)
            except KeyError as e:
                neg_dict[(cur_chrom, start)] = [read]
        else:
            pos_list.append(read)

        prev_pos = cur_pos
        prev_chrom = cur_chrom

    for x in neg_dict.keys():
        collapse_base(neg_dict[x], outBam, randomer, total_count, removed_count, max_hamming_distance, em)
        del neg_dict[x]

    collapse_base(pos_list, outBam, randomer, total_count, removed_count, max_hamming_distance, em)
    
    out_pos = prev_pos
    out_total.write("\t".join(map(str, [out_pos, cur_count])) + "\n")
    out_barcodes.write("\t".join(map(str, [out_pos, len(barcode_set)])) + "\n")

    inBam.close()
    outBam.close()
    out_total.close()
    out_barcodes.close()
    return total_count, removed_count


def output_metrics(metrics_file):
    with open(metrics_file, 'w') as metrics:
        metrics.write("\t".join(["randomer", "total_count", "removed_count"]) + "\n")
        for barcode in total_count.keys():
            metrics.write("\t".join(map(str, [barcode, total_count[barcode], removed_count[barcode]])) + "\n")


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bam", dest="bam", help="bam file to barcode collapse")
    parser.add_option("-c", "--randomer", action="store_true", dest="randomer", help="bam files are iclip barcoded")
    parser.add_option("-o", "--out_file", dest="out_file")
    parser.add_option("-m", "--metrics_file", dest="metrics_file")
    parser.add_option("-d", "--max_hamming_distance", dest="max_hamming_distance", default=0)
    parser.add_option("-e", "--em", action="store_true", default=False)

    (options, args) = parser.parse_args()

    
    if not (options.bam.endswith(".bam")):
        raise TypeError("%s, not bam file" % options.bam)

    pysam.index(options.bam)
    total_count, removed_count = barcode_collapse(options.bam, 
                                                  options.out_file,
                                                  options.randomer,
                                                  options.max_hamming_distance,
                                                  options.em
                                                  )

    output_metrics(options.metrics_file)
    pysam.index(options.out_file)
    sys.exit(0)
