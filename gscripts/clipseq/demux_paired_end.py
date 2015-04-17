__author__ = 'gpratt'


"""

Converts randomer + barcoded fastq files into something that can be barcode collapsed and mapped

"""

from collections import Counter
import gzip
import os
from optparse import OptionParser


def reformat_read(name_1, seq_1, plus_1, quality_1,
                  name_2, seq_2, plus_2, quality_2, barcodes,
                  RANDOMER_LENGTH=2):
    """ reformats read to have correct barcode attached
        name - read name
        seq - read sequence
        plus - +
        quality - read quality sequence this is a poor mans datastructure, designed for speed
        barcodes, dictionary of barcodes to search for

        returns str - barcode barcode found, str - randomer identified, str - reformateed read
    """
    barcode_lengths = {len(barcode) for barcode in barcodes.keys()}
    barcode = None

    #assumes larger barcodes are less likely, and searches for them first
    #for each barcode length see if known barcodes appear
    for barcode_length in sorted(barcode_lengths, reverse=True):
        cur_barcode = seq_1[:barcode_length]
        randomer = seq_2[:RANDOMER_LENGTH]
        if cur_barcode in barcodes:
            barcode = cur_barcode
            break

    name_1 = name_1[0] + randomer + ":" + name_1[1:]
    seq_1 = seq_1[barcode_length:]
    quality_1 = quality_1[barcode_length:]

    name_2 = name_2[0] + randomer + ":" + name_2[1:]
    seq_2 = seq_2[RANDOMER_LENGTH:]
    quality_2 = quality_2[RANDOMER_LENGTH:]

    #if none appear the barcode is unassigne
    if barcode is None:
        barcode = "unassigned"

    result_1 = name_1 + seq_1 + plus_1 + quality_1
    result_2 = name_2 + seq_2 + plus_2 + quality_2

    return barcode, randomer, result_1, result_2

if __name__ == "__main__":
    usage = """ takes raw fastq files and demultiplex inline randomer + adapter sequences  """
    parser = OptionParser(usage)
    parser.add_option("--fastq_1", dest="fastq_1", help="fastq file to barcode seperate")
    parser.add_option("--fastq_2", dest="fastq_2", help="fastq file to barcode seperate")

    parser.add_option("-b", "--barcodes", dest="barcodes", help="file of barcode / barcode id (tab sepearted, one barcode / barcode id on each line")
    parser.add_option("--out_file_1", dest="out_file_1")
    parser.add_option("--out_file_2", dest="out_file_2")
    parser.add_option("--length", type=int, dest="length", help="Number of randomers on the front of the second read", default=3)

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
            split_file_1.insert(-1, line[1])

            split_file_2 = options.out_file_2.split(".")
            split_file_2.insert(-1, line[1])

            barcodes[line[0]] = [gzip.open(".".join(split_file_1), 'w'),
                                 gzip.open(".".join(split_file_2), 'w'),
                                 ]

            randomer_counts[line[0]] = Counter()

    split_file_1 = options.out_file_1.split(".")
    split_file_1.insert(-1, "unassigned")

    split_file_2 = options.out_file_2.split(".")
    split_file_2.insert(-1, "unassigned")

    barcodes['unassigned'] = [gzip.open(".".join(split_file_1), 'w'),
                              gzip.open(".".join(split_file_2), 'w'),
                              ]
    randomer_counts['unassigned'] = Counter()

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

                barcode, randomer, result_1, result_2 = reformat_read(name_1,
                                                          seq_1,
                                                          plus,
                                                          quality_1,
                                                          name_2,
                                                          seq_2,
                                                          plus,
                                                          quality_2,
                                                          barcodes,
                                                          RANDOMER_LENGTH
                                                          )
                randomer_counts[barcode][randomer] += 1
                barcodes[barcode][0].write(result_1)
                barcodes[barcode][1].write(result_2)
            except StopIteration:
                break
        for barcode, randomers in randomer_counts.items():
            for randomer, count in randomers.items():
                metrics_file.write("%s\t%s\t%s\n" % (barcode, randomer, count))

    #cleans up at the end
    for lst in barcodes.values():
        for fn in lst:
            fn.close()

