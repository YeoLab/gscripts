__author__ = 'olga'

import pandas as pd

def max_csv(x):
    '''
    Of integers separated by commas, take the max
    e.g. 75112684,75112684 would return 75112684
    or 75112684,75112689 would return 75112689
    '''
    return max(map(int, x.split(',')))

def min_csv(x):
    '''
    Of integers separated by commas, take the minimum
    e.g. 75112684,75112684 would return 75112684
    or 75112684,75112689 would return 75112684
    '''
    return min(map(int, x.split(',')))


def read_miso_summary(filename):
    '''
    Reads a MISO summary file as a pandas dataframe, and adds a
    copy-paste-able genome location at the end, based on the minimum
    mRNA_starts and maximum mRNA_ends.

    '''
    df = pd.read_table(filename)
    genome_location = pd.DataFrame([ '%s:%d-%d' % (chrom, min_csv(starts), max_csv(stops))
                   for chrom, starts, stops in zip(df.chrom,
                                                   df.mRNA_starts,
                                                   df.mRNA_ends)],
                               columns=['genome_location'])
    return pd.concat([df, genome_location], axis=1)