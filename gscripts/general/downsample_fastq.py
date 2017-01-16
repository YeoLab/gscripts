__author__ = 'gpratt'
__author__ = 'gpratt'

__author__ = 'gpratt'
import argparse
import subprocess
import os


def wrap_wait_error(wait_result):
    if wait_result != 0:
        raise NameError("Failed to execute command correctly {}".format(wait_result))

def pre_process_fastq(fq, fq01, fq02, fq03, fq04, fq05, fq06, fq07, fq08, fq09):
    #split bam file into two, return file handle for the two bam files

    p = subprocess.Popen("zcat {} | wc -l".format(fq), shell=True, stdout=subprocess.PIPE) # Number of reads in the tagAlign file
    stdout, stderr = p.communicate()
    nlines = int(stdout.split()[0])
    print nlines
    fq_and_percent = [(fq01, int(nlines * .1)),
                      (fq02, int(nlines * .2)),
                      (fq03, int(nlines * .3)),
                      (fq04, int(nlines * .4)),
                      (fq05, int(nlines * .5)),
                      (fq06, int(nlines * .6)),
                      (fq07, int(nlines * .7)),
                      (fq08, int(nlines * .8)),
                      (fq09, int(nlines * .9))]

    cmds = []
    for fq_file, cur_lines in fq_and_percent:
        print cur_lines
        cur_lines -= cur_lines % 4 #Fixes issue where length might not be devisible by 4
        print "zcat {0} | head -n {1} | gzip > {2}".format(fq, cur_lines, fq_file)
        p1 = subprocess.Popen("zcat {0} | head -n {1} | gzip > {2}".format(fq, cur_lines, fq_file), shell=True)
        wrap_wait_error(p1.wait())

    return fq01, fq02

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Downsamples bam to a given number of reads')
    parser.add_argument(
        '--fastq', required=True, help='bam file to split')

    parser.add_argument(
        '--fq01', required=True, help='name of first output bam')
    parser.add_argument(
        '--fq02', required=True, help='name of second output bam')

    parser.add_argument(
        '--fq03', required=True, help='name of third output bam')

    parser.add_argument(
        '--fq04', required=True, help='name of fourth output bam')

    parser.add_argument(
        '--fq05', required=True, help='name of fifth output bam')

    parser.add_argument(
        '--fq06', required=True, help='name of sixth output bam')

    parser.add_argument(
        '--fq07', required=True, help='name of seventh output bam')

    parser.add_argument(
        '--fq08', required=True, help='name of eighth  output bam')

    parser.add_argument(
        '--fq09', required=True, help='name of ninth output bam')

    args = parser.parse_args()

    bam01, bam02 = pre_process_fastq(args.fastq, args.fq01, args.fq02, args.fq03,
                                   args.fq04, args.fq05, args.fq06, args.fq07,
                                   args.fq08, args.fq09)
