__author__ = 'gpratt'
import argparse

import subprocess


def sort_fastq_inplace(fastq, out_fastq):
    call = "cat {0} | sed  's/\t/ /' | paste  - - - - | sort -k1,1 -S 3G | tr '\t' '\n' > foo.fq && mv foo.fq {1}".format(fastq, out_fastq)
    subprocess.call(call, shell=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""Sorts in memory a single (small) fastq file.  This is important
    because STAR sometimes doesn't properly sort the output mate1 and mate2 files so we need a programatic control to keep this
    error from happening""")
    parser.add_argument("--in_fastq_1", help="Fastq File To Sort", required=True)
    parser.add_argument("--out_fastq_1", help="Sorted Fastq", required=True)
    parser.add_argument("--in_fastq_2", help="Fastq File To Sort", required=True)
    parser.add_argument("--out_fastq_2", help="Sorted Fastq", required=True)

    args = parser.parse_args()
    sort_fastq_inplace((args.inFastq_1, args.outFastq_1))
    sort_fastq_inplace((args.inFastq_2, args.outFastq_2))


#If this solution stops working revisit this post
#https://www.biostars.org/p/15011/