"""

Converts randomer + barcoded fastq files into something that can be barcode collapsed and mapped

"""

from collections import Counter
import gzip
import os
from optparse import OptionParser


def reformat_read(name, seq, plus, quality, barcodes,
                  RANDOMER_FRONT_LENGTH=3, RANDOMER_BACK_LENGTH=2):
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
        adapter_length = RANDOMER_FRONT_LENGTH + RANDOMER_BACK_LENGTH + barcode_length
        randomer = seq[:adapter_length]
        cur_barcode = randomer[RANDOMER_FRONT_LENGTH:adapter_length - RANDOMER_BACK_LENGTH]

        if cur_barcode in barcodes:
            barcode = cur_barcode
            break
    
    name = name[0] + randomer + ":" + name[1:]
    seq = seq[adapter_length:]
    quality = quality[adapter_length:]
    
    #if none appear the barcode is unassigne
    if barcode is None:
        barcode = "unassigned"
        
    result = name + seq + plus + quality
    return barcode, randomer, result

if __name__ == "__main__":
    usage = """ takes raw fastq files and demultiplex inline randomer + adapter sequences  """
    parser = OptionParser(usage)
    parser.add_option("-f", "--fastq", dest="fastq", help="fastq file to barcode seperate")
    parser.add_option("-b", "--barcodes", dest="barcodes", help="file of barcode / barcode id (tab sepearted, one barcode / barcode id on each line")
    parser.add_option("-o", "--out_file", dest="out_file")
    parser.add_option("-m", "--metrics_file", dest="metrics_file")
    parser.add_option("--front", type=int, dest="front_length", help="Number of randomers before the barcode", default=3)
    parser.add_option("--back", type=int, dest="back_length", help="Number of randomers after the barcode", default=2)
    
    (options,args) = parser.parse_args()
    
    #if a gziped file then we reassign open to be the gzip open and continue using that for the rest of the
    #program
    my_open = gzip.open if os.path.splitext(options.fastq)[1] == ".gz" else open
    #creates different barcode files to assign reads to

    RANDOMER_FRONT_LENGTH = options.front_length
    RANDOMER_BACK_LENGTH = options.back_length

    barcodes = {}
    randomer_counts = {} 
    with open(options.barcodes) as barcodes_file:
        for line in barcodes_file:
            line = line.strip("\n").split("\t")
            split_file = options.out_file.split(".")
            split_file.insert(-2, line[1])
            barcodes[line[0]] = gzip.open(".".join(split_file), 'w')
            randomer_counts[line[0]] = Counter()
    
    split_file = options.out_file.split(".")
    split_file.insert(-2, "unassigned")
    barcodes['unassigned'] = gzip.open(".".join(split_file), 'w')
    randomer_counts['unassigned'] = Counter()
    
    #reads through initial file parses everything out
    with my_open(options.fastq) as fastq_file, open(options.metrics_file, 'w') as metrics_file:
        while True:
            try:
                name = fastq_file.next()
                seq = fastq_file.next()
                fastq_file.next() #got to consume the read
                plus = "+\n" #sometimes the descriptor is here, don't want it
                quality = fastq_file.next()
                
                barcode, randomer, result = reformat_read(name, seq, plus, quality, barcodes,
                                                          RANDOMER_FRONT_LENGTH, RANDOMER_BACK_LENGTH)
                randomer_counts[barcode][randomer] += 1
                barcodes[barcode].write(result)
            except StopIteration:
                break
        for barcode, randomers in randomer_counts.items():
            for randomer, count in randomers.items():
                metrics_file.write("%s\t%s\t%s\n" % (barcode, randomer, count))
    
    #cleans up at the end
    for fn in barcodes.values():
        fn.close()

