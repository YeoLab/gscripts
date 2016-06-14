__author__ = 'gpratt'

import sys
import argparse
from collections import Counter

import pyBigWig
import numpy as np
from scipy import stats
import pybedtools


def big_wig_corr(full, semi, regions):

    full = pyBigWig.open(full)
    semi = pyBigWig.open(semi)
    regions = pybedtools.BedTool(regions)

    full_result = []
    semi_result = []
    for interval in regions:
        gene_full_values = np.array(full.values(interval.chrom, interval.start, interval.stop))
        gene_semi_values = np.array(semi.values(interval.chrom, interval.start, interval.stop))

        filtered_gene_full_values = gene_full_values[~np.isnan(gene_full_values) & (gene_full_values != 0)]
        filtered_gene_semi_values = gene_semi_values[~np.isnan(gene_full_values) & (gene_full_values != 0)]
        filtered_gene_semi_values = np.nan_to_num(filtered_gene_semi_values)

        full_result.append(filtered_gene_full_values)
        semi_result.append(filtered_gene_semi_values)
    full_result = np.concatenate(full_result)
    semi_result = np.concatenate(semi_result)

    return stats.pearsonr(full_result, semi_result)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Correlates two bigwig files with pearosnR this is different from other tools because it allows for '
                    'zeros in one of the files (the semi file).  The full file does not allow for zeros still')
    parser.add_argument('--full', help='input full bw file')
    parser.add_argument('--semi', help='input semi bw file')
    parser.add_argument('--regions', help='input bed file of regions to run correlation on')

    args = parser.parse_args()

    results = big_wig_corr(args.full, args.semi, args.regions)
    print args.full, args.semi, results[0]

    sys.exit(0)