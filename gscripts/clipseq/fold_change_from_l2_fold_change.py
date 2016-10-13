__author__ = 'gpratt'

import argparse

import pybedtools
import sys

def fix_score(interval):
    interval.score = str(pow(2, float(interval.score)))
    return interval

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="fixed log2 fold change enrichment to just be fold enrichment for a compressed input norm peak file")

    parser.add_argument(
        '--in_file', required=True, help='name of  compressed input normed peak file')

    parser.add_argument(
        '--out_file', required=True, help='name output file')

    args = parser.parse_args()

    rep1_out = pybedtools.BedTool(args.in_file).each(fix_score).saveas(args.out_file)
    sys.exit(0)

