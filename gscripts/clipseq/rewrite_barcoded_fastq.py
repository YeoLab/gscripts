from collections import Counter
import gzip
from optparse import OptionParser
import sys

from Bio import SeqIO

parser = OptionParser()
parser.add_option("-f", "--fastq", dest="fastq", help="fastq file to barcode seperate")
parser.add_option("-b", "--barcodes", dest="barcodes", help="file of barcode / barcode id")
parser.add_option("-o", "--out_file", dest="out_file")
parser.add_option("-m", "--metrics_file", dest="metrics_file")

(options,args) = parser.parse_args()

#reformats read, writes to correct barcode file if exists, else writes to unassigned file
#records barcode stats
def reformat_read(file_handle, barcodes):
    name = file_handle.next()
    seq = file_handle.next()
    plus = file_handle.next()
    quality = file_handle.next()
    
    randomer = seq[:9]
    name = name[0] + randomer + ":" + name[1:]
    seq = seq[9:]
    quality = quality[9:]
    barcode = randomer[3:7]

    if barcode not in barcodes:
        barcode = "unassigned"
        
    randomer_counts[barcode][randomer] += 1

    barcodes[barcode].write(name)
    barcodes[barcode].write(seq)
    barcodes[barcode].write(plus)
    barcodes[barcode].write(quality)
    

#creates different barcode files to assign reads to
barcodes = {}
randomer_counts = {} 
with open(options.barcodes) as barcodes_file:
    for line in barcodes_file:
        line = line.strip().split()
        print options.out_file
        split_file = options.out_file.split(".")
        split_file.insert(-1, line[1])
        barcodes[line[0]] = open(".".join(split_file), 'w')
        randomer_counts[line[0]] = Counter()

split_file = options.out_file.split(".")
split_file.insert(-1, "unassigned")
barcodes['unassigned'] = open(".".join(split_file), 'w')
randomer_counts['unassigned'] = Counter()

#reads through initial file parses everything out
with open(options.fastq) as fastq_file, open(options.metrics_file, 'w') as metrics_file:
    while True:
        try:
            reformat_read(fastq_file, barcodes)
        except StopIteration:
            break
    for barcode, randomers in randomer_counts.items():
        for randomer, count in randomers.items():
            metrics_file.write("%s\t%s\t%s\n" % (barcode, randomer, count))

#cleans up at the end
for file in barcodes.values():
    file.close()

