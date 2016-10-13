from collections import Counter, defaultdict, OrderedDict
from itertools import izip
import gzip
import os
from optparse import OptionParser

RANDOMER_LENGTH = 10
if __name__ == "__main__":
    usage = """ takes raw fastq files and demultiplex inline randomer + adapter sequences  """
    parser = OptionParser(usage)
    parser.add_option("--fastq1", dest="fastq1", help="fastq file to barcode seperate")
    parser.add_option("--fastq2", dest="fastq2", help="fastq file to barcode seperate")
    parser.add_option("--out_file1", dest="out_file1")
    parser.add_option("--out_file2", dest="out_file2")


    (options, args) = parser.parse_args()

    my_open = gzip.open if os.path.splitext(options.fastq1)[1] == ".gz" else open

    out_file1 = my_open(options.out_file1, 'w')
    out_file2 = my_open(options.out_file2, 'w')

    with my_open(options.fastq1) as fastq_file1, my_open(options.fastq2) as fastq_file2:
        while True:
            try:
                plus = "+\n" #sometimes the descriptor is here, don't want it

                name_1 = fastq_file1.next()
                seq_1 = fastq_file1.next()
                fastq_file1.next() #got to consume the read
                quality_1 = fastq_file1.next()


                name_2 = fastq_file2.next()
                seq_2 = fastq_file2.next()
                fastq_file2.next() #got to consume the read
                quality_2 = fastq_file2.next()

                randomer = seq_1[:RANDOMER_LENGTH]

                name_1 = name_1[0] + randomer + ":" + name_1[1:]
                seq_1 = seq_1[RANDOMER_LENGTH:]
                quality_1 = quality_1[RANDOMER_LENGTH:]

                name_2 = name_2[0] + randomer + ":" + name_2[1:]

                out_file1.write(name_1 + seq_1 + "+\n" + quality_1)
                out_file2.write(name_2 + seq_2 + "+\n" + quality_2)

            except StopIteration:
                break
