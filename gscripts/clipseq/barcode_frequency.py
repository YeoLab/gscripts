'''
Created on Dec 11, 2013

@author: gpratt
'''

from collections import defaultdict
import gzip
from optparse import OptionParser
import os

import numpy as np
import pandas as pd

def handle_seq(seq, barcode_df, result_dict):
    """
    
    seq - dna sequence to look for barcodes in
    barcode_df - pandas dataframe of barcodes, ids
    result_dict - dictionary where results are stored, not returned for speed
    
    counts number of barcodes in each location for each read
    
    """
    for i in range(len(seq)):
        for barcode in barcode_df.index:
            possible_match = seq[i: i + len(barcode)]
            if possible_match == barcode:
                #old style, default dicts are slow, could make into matrix
                if i not in result_dict[barcode]: 
                    result_dict[barcode][i] = 0
                result_dict[barcode][i] += 1

def calculate_barcode_frequency(file_handle, barcodes):
    
    """
    
    fastq -- file handle to fastq test text file (for fast processing)
    barcodes -- csv file handle (or file) of barcode\tbarcode_id
    
    returns pandas dataframe of frequency of each barcode at each base
    
    """
    
    barcodes_df = pd.read_csv(barcodes, 
                              sep="\t", 
                              names= ["seq", "id"], 
                              index_col=0)
    
    result_dict = defaultdict(dict)
    
    while True:
        file_handle.next() #name
        seq = file_handle.next()
        file_handle.next() #plus
        file_handle.next() #qual
        handle_seq(seq, barcodes, result_dict)
        
    data_frame = pd.DataFrame(result_dict).T
    data_frame[np.isnan(data_frame)] = 0
    data_frame.index.name = "barcode"
    data_frame.insert(0, "barcode_name", barcodes_df.id)
    
    return data_frame

if __name__ == '__main__':
    usage = """ calculates frequency of each barcode sequence across all reads for
    all barcodes """
    
    parser = OptionParser(usage)
    parser.add_option("-f", "--fastq", dest="fastq", help="fastq file to barcode seperate")
    parser.add_option("-b", "--barcodes", dest="barcodes", help="file of barcode / barcode id")
    parser.add_option("-o", "--out_file", dest="out_file")
    
    (options,args) = parser.parse_args()
    if os.path.splitext(options.fastq)[1] == ".gz":
        handle = gzip.open(options.fastq)
    else:
        handle = open(options.fastq)
        
    calculate_barcode_frequency(handle, options.barcodes).to_csv(options.out_file)
    
 