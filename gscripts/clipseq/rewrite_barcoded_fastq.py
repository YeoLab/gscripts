from collections import Counter
import gzip
from optparse import OptionParser
import sys

from Bio import SeqIO

parser = OptionParser()
parser.add_option("-f", "--fastq", dest="fastq", help="fastq file to barcode seperate")
parser.add_option("-b", "--barcode", dest="barcode", default=.05, help="barcode sequence to split on")
parser.add_option("-o", "--out_file", dest="out_file")
parser.add_option("-m", "--metrics_file", dest="metrics_file")

(options,args) = parser.parse_args()
barcode_counts = Counter()

def reformat_read(read):
    read.id = str(read.seq[:9]) + ":" + read.name
    barcode_counts[str(read.seq[:9])] += 1
    
    return read[9:]

with open(options.fastq) as fastq_file, open(options.out_file, 'w') as out_file, open(options.metrics_file, 'w') as metrics_file:
    read_iterator = (reformat_read(read) for read in SeqIO.parse(fastq_file, 'fastq'))

    SeqIO.write(read_iterator, out_file, 'fastq')
    for barcode, count in sorted(barcode_counts.items(), key=lambda x: x[1], reverse=True):
        metrics_file.write("%s\t%s\n" % (barcode, count))
