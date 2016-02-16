__author__ = 'gpratt'


"""

Converts randomer + barcoded fastq files into something that can be barcode collapsed and mapped

"""

from collections import Counter, defaultdict
from itertools import izip
import gzip
import os
from optparse import OptionParser


def hamming(word1, word2):
    """
    Gets hamming distance between two words, this is an odd implementation because
    we actually want Ns to appear similar to other barcodes as this makes the results more stringent, not less
    :param word1:
    :param word2:
    :return:
    """
    return sum(a != b and not (a == "N" or b == "N") for a, b in izip(word1, word2))


def read_has_barcode(barcodes, read, max_hamming_distance=0):
    """
    Checks if read has one of the barcodes listed in barcodes list and returns back that barcode

    :param barcodes: list of barcode sequence
    :param read: str, read to check if barcode exists in
    :param hamming: max hamming distance between given barcode and barcode in read
    :return: returns best barcode match in read, none if none found
    """

    #This takes care of the edge case of some barcodes having Ns on the end
    effective_barcodes = {barcode.split("N")[0]: barcode for barcode in barcodes}

    barcode_lengths = {len(barcode) for barcode in effective_barcodes}
    closest_barcode = None
    #assumes larger barcodes are less likely, and searches for them first
    #for each barcode length see if known barcodes appear
    for barcode_length in sorted(barcode_lengths, reverse=True):
        #Gets min hamming distance between barcode in read and barcode in list, behavior undefined if tie
        read_barcode = read[:barcode_length]
        cur_length_barcodes = [barcode for barcode in effective_barcodes if len(barcode) == barcode_length]
        hamming_distances = {barcode: hamming(barcode, read_barcode) for barcode in cur_length_barcodes}
        min_barcode = min(hamming_distances, key=hamming_distances.get)
        if hamming_distances[min_barcode] <= max_hamming_distance:
            closest_barcode = effective_barcodes[min_barcode]
            break
    return closest_barcode, read_barcode


def reformat_read(name_1, seq_1, plus_1, quality_1,
                  name_2, seq_2, plus_2, quality_2, barcodes,
                  RANDOMER_LENGTH=2, max_hamming_distance=0):
    """ reformats read to have correct barcode attached
        name - read name
        seq - read sequence
        plus - +
        quality - read quality sequence this is a poor mans datastructure, designed for speed
        barcodes, dictionary of barcodes to search for

        returns str - barcode barcode found, str - randomer identified, str - reformateed read
    """


    barcode, actual_barcode = read_has_barcode(barcodes.keys(), seq_1, max_hamming_distance)
    barcode_length = len(barcode) if barcode is not None else 0

    randomer = seq_2[:RANDOMER_LENGTH]

    name_1 = name_1[0] + randomer + ":" + name_1[1:]
    seq_1 = seq_1[barcode_length:]
    quality_1 = quality_1[barcode_length:]

    name_2 = name_2[0] + randomer + ":" + name_2[1:]
    seq_2 = seq_2[RANDOMER_LENGTH:]
    quality_2 = quality_2[RANDOMER_LENGTH:]

    #if none appear the barcode is unassigned
    if barcode is None:
        barcode = "unassigned"

    result_1 = name_1 + seq_1 + plus_1 + quality_1
    result_2 = name_2 + seq_2 + plus_2 + quality_2

    return barcode, actual_barcode, randomer, result_1, result_2

if __name__ == "__main__":
    usage = """ takes raw fastq files and demultiplex inline randomer + adapter sequences  """
    parser = OptionParser(usage)
    parser.add_option("--fastq_1", dest="fastq_1", help="fastq file to barcode seperate")
    parser.add_option("--fastq_2", dest="fastq_2", help="fastq file to barcode seperate")

    parser.add_option("-b", "--barcodes", dest="barcodes", help="file of barcode / barcode id (tab sepearted, one barcode / barcode id on each line")
    parser.add_option("--out_file_1", dest="out_file_1")
    parser.add_option("--out_file_2", dest="out_file_2")
    parser.add_option("--length", type=int, dest="length", help="Number of randomers on the front of the second read", default=3)
    parser.add_option("--max_hamming_distance", type=int, dest="max_hamming_distance", help="Max Hamming distance between read barcode and given barcodes to assign a read to a given barcode", default=1)

    parser.add_option("-m", "--metrics_file", dest="metrics_file")

    (options,args) = parser.parse_args()

    #if a gziped file then we reassign open to be the gzip open and continue using that for the rest of the
    #program
    my_open = gzip.open if os.path.splitext(options.fastq_1)[1] == ".gz" else open
    #creates different barcode files to assign reads to

    RANDOMER_LENGTH = options.length

    barcodes = {}
    randomer_counts = {}
    with open(options.barcodes) as barcodes_file:
        for line in barcodes_file:
            line = line.strip().split("\t")
            split_file_1 = options.out_file_1.split(".")
            split_file_1.insert(-2, line[1])

            split_file_2 = options.out_file_2.split(".")
            split_file_2.insert(-2, line[1])

            barcodes[line[0]] = [gzip.open(".".join(split_file_1), 'w'),
                                 gzip.open(".".join(split_file_2), 'w'),
                                 ]

            randomer_counts[line[0]] = defaultdict(Counter)

    split_file_1 = options.out_file_1.split(".")
    split_file_1.insert(-2, "unassigned")

    split_file_2 = options.out_file_2.split(".")
    split_file_2.insert(-2, "unassigned")

    barcodes['unassigned'] = [gzip.open(".".join(split_file_1), 'w'),
                              gzip.open(".".join(split_file_2), 'w'),
                              ]
    randomer_counts['unassigned'] = defaultdict(Counter)

    #reads through initial file parses everything out
    with my_open(options.fastq_1) as fastq_file_1, my_open(options.fastq_2) as fastq_file_2, open(options.metrics_file, 'w') as metrics_file:
        while True:
            try:
                name_1 = fastq_file_1.next()
                seq_1 = fastq_file_1.next()
                fastq_file_1.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality_1 = fastq_file_1.next()

                name_2 = fastq_file_2.next()
                seq_2 = fastq_file_2.next()
                fastq_file_2.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality_2 = fastq_file_2.next()
                if name_1.split()[0] != name_2.split()[0]:
                    print name_1, name_2
                    raise Exception("Read 1 is not same name as Read 2")

                barcode, actual_barcode, randomer, result_1, result_2 = reformat_read(name_1, seq_1, plus, quality_1,
                                                                      name_2, seq_2, plus, quality_2,
                                                                      barcodes, RANDOMER_LENGTH,
                                                                      max_hamming_distance=options.max_hamming_distance)
                randomer_counts[barcode][actual_barcode][randomer] += 1
                barcodes[barcode][0].write(result_1)
                barcodes[barcode][1].write(result_2)
            except StopIteration:
                break
        for barcode, actual_barcodes in randomer_counts.items():
            for actual_barcode, randomers in actual_barcodes.items():
                for randomer, count in randomers.items():
                    metrics_file.write("%s\t%s\t%s\t%s\n" % (barcode, actual_barcode, randomer, count))

    #cleans up at the end
    for lst in barcodes.values():
        for fn in lst:
            fn.close()

