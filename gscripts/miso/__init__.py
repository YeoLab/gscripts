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
    Reads a MISO summary file as a pandas dataframe, and adds these columns:
    1. a copy-paste-able genome location at the end, based on the minimum
       mRNA_starts and maximum mRNA_ends. (df.genome_location)
    2. The difference between df.ci_high and df.ci_low (df.ci_diff)
    3. The left and right halves of the confidence interval, e.g. the right
       half is df.ci_high - df.miso_posterior_mean. (df.ci_left_half and
       df.ci_right_half)
    4. The max of the two left and right confidence interval halves
       (df.ci_halves_max)
    '''
    df = pd.read_table(filename, index_col=0)
    genome_location = pd.DataFrame(
        ['%s:%d-%d' % (chrom, min_csv(starts), max_csv(stops))
         for chrom, starts, stops in zip(df.chrom,
                                         df.mRNA_starts,
                                         df.mRNA_ends)],
        columns=['genome_location'], index=df.index)
    ci_diff = pd.DataFrame(df.ci_high - df.ci_low, columns=['ci_diff'],
                           index=df.index)
    ci_halves = pd.DataFrame(
        {'ci_left_half': (df.ci_high - df.miso_posterior_mean),
         'ci_right_half': (df.miso_posterior_mean - df.ci_low)},
        index=df.index)
    ci_halves_max = pd.DataFrame(ci_halves.max(axis=1),
                                        columns=['ci_halves_max'])
    return pd.concat([df, genome_location, ci_diff, ci_halves,
                      ci_halves_max], axis=1)