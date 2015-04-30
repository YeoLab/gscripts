#!/usr/bin/env python

__author__ = 'lovci'

import argparse
import subprocess
import os


def quality_processing_call(fastq):
    return " ".join(['cutadapt', '-q', '37', '-f', 'fastq', '\'' + fastq + '\''])


def construct_cutadapt_call(miRNA_seq, miRNA_id, output, **kwargs):
    miRNA_seq = miRNA_seq.replace("U", "T")
    return " ".join(['cutadapt',
                     '-g', 'N' + miRNA_seq,
                     '-g', 'N' + miRNA_seq[1:],
                     '-g', 'N' + miRNA_seq[2:],
                     '-g', 'N' + miRNA_seq[3:],
                     '-g', 'N' + miRNA_seq[1:-1],
                     '-g', 'N' + miRNA_seq[1:-2],
                     '-g', 'N' + miRNA_seq[:-2],
                     '--match-read-wildcards', '--trimmed-only',
                     '-m', '18', '-y', '\'.' + miRNA_id + '\'', '-f', 'fastq', '-O', '20', '-e', '0.05',
                     '-o', '\'' + output + '\'', '-'])


def run_cutadapt(miRNA_seq, miRNA_id, fastq, output, **kwargs):
    call_quality = quality_processing_call(fastq)
    call_cutadapt = construct_cutadapt_call(miRNA_seq, miRNA_id, output)
    server1 = subprocess.Popen(call_quality, shell=True, stdout=subprocess.PIPE)
    server2 = subprocess.Popen(call_cutadapt, shell=True, stdin=server1.stdout)
    server1.wait()
    server2.wait()


def count_lines(file_name):
    line_count = int(subprocess.check_output(['wc', '-l', file_name]).split(' ')[0])
    return line_count


parser = argparse.ArgumentParser(description='split a fastq file')
parser.add_argument('--miRNA_seq', dest='miRNA_seq', type=str,
                   help='DNA or RNA sequence of miRNA')

parser.add_argument('--miRNA_id', dest='miRNA_id', type=str,
                   help='name of miRNA')

parser.add_argument('--fastq', dest='fastq', type=str,
                   help='name of fastq to split')

parser.add_argument('--output', dest='output', type=str,
                   help='name of output fastq, after split')

parser.add_argument('--print_call', dest='print_call',
                    action='store_true', default=False,
                    help='only prints the call, doesn\'t actually run cutadapt')

args = vars(parser.parse_args())

if args['print_call']:
    print construct_cutadapt_call(**args)
else:
    run_cutadapt(**args)
    lc = count_lines(args['output'])

    try:
        assert lc >= 4 #one fastq entry
        print "output for %s on %s is %d lines long" %(args['miRNA_id'], args['fastq'], lc/4)
    except:
        os.remove(args['output'])
        raise AssertionError("output was too small for %s on %s, raising a failure status" % (args['miRNA_id'],
                                                                                            args['fastq']))
