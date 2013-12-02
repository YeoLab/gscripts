""""

fixes scores so we can actually make a bigbed track

"""

from optparse import OptionParser


import pybedtools


def adjust_score(read):
    import numpy as np
    if float(read.score) != 0.0:
        read.score = str(int(-10 * np.log10(float(read.score))))
        if int(read.score) > 1000:
            read.score = str(1000)
    else:
        read.score = str(1000)
    return read
        
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-b", "--bed", dest="bed", help="bed file to barcode collapse")
    parser.add_option("-o", "--out_file", dest="out_file")
    
    
    
    (options,args) = parser.parse_args()
    
    pybedtools.BedTool(options.bed).each(adjust_score).saveas(options.out_file)
    

