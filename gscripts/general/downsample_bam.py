__author__ = 'gpratt'

__author__ = 'gpratt'
import argparse
import subprocess
import os


def wrap_wait_error(wait_result):
    if wait_result != 0:
        raise NameError("Failed to execute command correctly {}".format(wait_result))

def pre_process_bam(bam, bam01, bam02, bam03, bam04, bam05, bam06, bam07, bam08, bam09, no_shuffle, no_sort):
    #split bam file into two, return file handle for the two bam files
    
    p = subprocess.Popen("samtools view {} | wc -l".format(bam), shell=True, stdout=subprocess.PIPE) # Number of reads in the tagAlign file
    stdout, stderr = p.communicate()
    nlines = int(stdout)
    
    p = subprocess.Popen("samtools view -H {} | wc -l".format(bam), shell=True, stdout=subprocess.PIPE) # Number of header lines (for when we've got a lot of chromosomes)
    stdout, stderr = p.communicate()
    n_header_lines = int(stdout)

    if no_shuffle:
        shuffled_bam = os.path.splitext(bam)[0]
    else: #shuffle
        shuffled_bam = os.path.splitext(os.path.basename(bam))[0] + "_shuff"    
        p = subprocess.Popen("samtools bamshuf {0} {1}".format(bam, shuffled_bam), shell=True) # This will shuffle the lines in the file and split it into two parts
        wrap_wait_error(p.wait())
    
    bam_and_percent = [(bam01, int(nlines * .1) + n_header_lines),
                       (bam02, int(nlines * .2) + n_header_lines),
                       (bam03, int(nlines * .3) + n_header_lines),
                       (bam04, int(nlines * .4) + n_header_lines),
                       (bam05, int(nlines * .5) + n_header_lines),
                       (bam06, int(nlines * .6) + n_header_lines),
                       (bam07, int(nlines * .7) + n_header_lines),
                       (bam08, int(nlines * .8) + n_header_lines), 
                       (bam09, int(nlines * .9) + n_header_lines),]

    cmds = []
    for bam_file, percent in bam_and_percent:
        if percent % 2 == 1:
            percent -= 1
        if no_sort: #if we are aren't shuffling, don't delete
        #Make sure I select pairs of reads
            cmd = "samtools view -h {0}.bam | head -n {1} | samtools view -bS - > {2}.bam".format(shuffled_bam, percent, os.path.splitext(bam_file)[0])
            p1 = subprocess.Popen(cmd, shell=True)
        else: #sort
            p1 = subprocess.Popen("samtools view -h {0}.bam | head -n {1} | samtools view -bS - | samtools sort - {2}  && samtools index {2}.bam".format(shuffled_bam, percent, os.path.splitext(bam_file)[0]), shell=True)
        wrap_wait_error(p1.wait())

    if not no_shuffle: #if we are aren't shuffling, don't delete
        p1 = subprocess.Popen("rm {0}.bam".format(shuffled_bam), shell=True)
        wrap_wait_error(p1.wait())

    return bam01, bam02

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Downsamples bam to a given number of reads')
    parser.add_argument(
        '--bam', required=True, help='bam file to split')


    parser.add_argument(
        '--bam01', required=True, help='name of first output bam')
    parser.add_argument(
        '--bam02', required=True, help='name of second output bam')

    parser.add_argument(
        '--bam03', required=True, help='name of third output bam')

    parser.add_argument(
        '--bam04', required=True, help='name of fourth output bam')

    parser.add_argument(
        '--bam05', required=True, help='name of fith output bam')

    parser.add_argument(
        '--bam06', required=True, help='name of sixth output bam')

    parser.add_argument(
        '--bam07', required=True, help='name of seventh output bam')

    parser.add_argument(
        '--bam08', required=True, help='name of eighth  output bam')

    parser.add_argument(
        '--bam09', required=True, help='name of ninth output bam')


    parser.add_argument("--no_shuffle", action="store_true", help="Don't shuffle input bam file, only use this if input bam is already somehow shuffled")
    parser.add_argument("--no_sort", action="store_true", help="Don't sort the resulting bam files")

    args = parser.parse_args()

    bam01, bam02 = pre_process_bam(args.bam, args.bam01, args.bam02, args.bam03,
                                   args.bam04, args.bam05, args.bam06, args.bam07,
                                   args.bam08, args.bam09, args.no_shuffle, args.no_sort)
