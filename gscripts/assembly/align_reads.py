__author__ = 'lovci'

trinity_home = "/home/mlovci/software/trinityrnaseq_r20131110"
read_align_util = os.path.join(trinity_home, "/util/alignReads.pl")

import argparse
parser = argparse.ArgumentParser(description='align reads with trinity\'s wrapper for gsnap')
parser.add_argument('--fastq', action='append')
parser.add_argument('--paired', default=False, action='store_true')

args = parser.parse_args()

if args.paired:
    #paired-end aligner
    pass
else:
    #single-end aligner