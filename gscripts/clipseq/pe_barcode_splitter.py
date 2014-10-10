__author__ = 'gpratt'

from collections import Counter
import gzip
import os
from optparse import OptionParser


def parse_barcode_file(barcodes_file, out_file1, out_file2):


    barcodes = {}
    with open(barcodes_file) as barcodes_file:
        for line in barcodes_file:
            line = line.strip("\n").split("\t")
            split_file1 = out_file1.split(".")
            split_file1.insert(-1, line[1])

            split_file2 = out_file2.split(".")
            split_file2.insert(-1, line[1])
            barcodes[line[0]] = [open(".".join(split_file1), 'w'), open(".".join(split_file2), 'w')]

    split_file1 = out_file1.split(".")
    split_file1.insert(-1, "unassigned")

    split_file2 = out_file2.split(".")
    split_file2.insert(-1, "unassigned")
    barcodes['unassigned'] = [open(".".join(split_file1), 'w'), open(".".join(split_file2), 'w')]

    return barcodes

def reformat_read(name, seq, plus, quality, barcodes):
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
        cur_barcode = seq[:barcode_length]
        if cur_barcode in barcodes:
            barcode = cur_barcode
            break

    seq = seq[barcode_length:]
    quality = quality[barcode_length:]

    #if none appear the barcode is unassigned
    if barcode is None:
        barcode = "unassigned"

    result = name + seq + plus + quality
    return barcode, result

if __name__ == "__main__":
    usage = """ takes raw fastq files and demultiplex inline randomer + adapter sequences  """
    parser = OptionParser(usage)
    parser.add_option("--fastq1", dest="fastq1", help="fastq file to barcode seperate")
    parser.add_option("--fastq2", dest="fastq2", help="fastq file to barcode seperate")
    parser.add_option("--barcodes", dest="barcodes", help="file of barcode / barcode id (tab sepearted, one barcode / barcode id on each line")
    parser.add_option("--out_file1", dest="out_file1")
    parser.add_option("--out_file2", dest="out_file2")

    parser.add_option("-m", "--metrics_file", dest="metrics_file")

    (options, args) = parser.parse_args()

    barcodes = parse_barcode_file(options.barcodes, options.out_file1, options.out_file2)

    my_open = gzip.open if os.path.splitext(options.fastq1)[1] == ".gz" else open

    with my_open(options.fastq1) as fastq_file1, my_open(options.fastq2) as fastq_file2, open(options.metrics_file, 'w') as metrics_file:
        while True:
            try:
                plus = "+\n" #sometimes the descriptor is here, don't want it

                name1 = fastq_file1.next()
                seq1 = fastq_file1.next()
                fastq_file1.next() #got to consume the read
                quality1 = fastq_file1.next()

                name2 = fastq_file2.next()
                seq2 = fastq_file2.next()
                fastq_file2.next() #got to consume the read
                quality2 = fastq_file2.next()

                barcode, result = reformat_read(name1, seq1, plus, quality1, barcodes)
                barcodes[barcode][0].write(result)
                barcodes[barcode][1].write(name2 + seq2 + plus + quality2)

            except StopIteration:
                break