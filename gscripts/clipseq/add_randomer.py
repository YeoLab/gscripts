__author__ = 'gpratt'

import gzip
from optparse import OptionParser
import sys


def add_back_randomers(in_file_name, out_file_name):
    #reads through initial file parses everything out
    with gzip.open(in_file_name) as fastq_file, gzip.open(out_file_name, 'w') as out_file:
        while True:
            try:
                name_1 = fastq_file.next()
                seq_1 = fastq_file.next()
                fastq_file.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality_1 = fastq_file.next()

                randomer = name_1.split(":")[0][1:]
                name_1 = "@" + ":".join(name_1.split(":")[1:])
                seq_1 = randomer + seq_1
                quality_1 = len(randomer) * "-" + quality_1

                out_file.write(name_1)
                out_file.write(seq_1)
                out_file.write(plus)
                out_file.write(quality_1)

            except StopIteration:
                break

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-i", "--in_file", help="bam file to add barcodes back to")
    parser.add_option("-o", "--out_file", help="out bam file")
    (options, args) = parser.parse_args()

    add_back_randomers(options.in_file, options.out_file)

    sys.exit(0)
