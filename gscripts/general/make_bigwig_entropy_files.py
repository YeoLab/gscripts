__author__ = 'gpratt'

import argparse
import subprocess
import os




if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Makes bigwig entropy tracks for an ip input pair')
    parser.add_argument(
        '--ip_bam', required=True, help='bam file to split')
    parser.add_argument(
        '--input_bam', required=True, help='name of first output bam')
    parser.add_argument("--bw_pos", help="positive bw file name", required=True)
    parser.add_argument("--bw_neg", help="negative bw file name", required=True)

    parser.add_argument(
        '--genome', required=True, help='name of second output bam')

    args = parser.parse_args()

    bam01, bam02 = pre_process_fastq(args.fastq, args.fq01, args.fq02, args.fq03,
                                   args.fq04, args.fq05, args.fq06, args.fq07,
                                   args.fq08, args.fq09)
