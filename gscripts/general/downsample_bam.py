__author__ = 'gpratt'

__author__ = 'gpratt'
import argparse
import subprocess
import os


def pre_process_bam(bam, bam01, bam02, bam03, bam04, bam05, bam06, bam07, bam08, bam09):
    #split bam file into two, return file handle for the two bam files
    p = subprocess.Popen("samtools view {} | wc -l".format(bam), shell=True, stdout=subprocess.PIPE) # Number of reads in the tagAlign file
    stdout, stderr = p.communicate()
    nlines = int(stdout)
    p = subprocess.Popen("samtools view {0} | shuf > {0}.shuffled.sam".format(bam), shell=True) # This will shuffle the lines in the file and split it into two parts
    p.wait()

    p1 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam01, int(nlines * .1)), shell=True)
    p2 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam02, int(nlines * .2)), shell=True)
    p3 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam03, int(nlines * .3)), shell=True)
    p4 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam04, int(nlines * .4)), shell=True)
    p5 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam05, int(nlines * .5)), shell=True)
    p6 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam06, int(nlines * .6)), shell=True)
    p7 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam07, int(nlines * .7)), shell=True)
    p8 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam08, int(nlines * .8)), shell=True)
    p9 = subprocess.Popen("head -n {2} {0}.shuffled.sam > {1}.tmp".format(bam, bam09, int(nlines * .9)), shell=True)

    p1.wait()
    p2.wait()
    p3.wait()
    p4.wait()
    p5.wait()
    p6.wait()
    p7.wait()
    p8.wait()
    p9.wait()

    p1 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam01), shell=True)
    p2 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam02), shell=True)
    p3 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam03), shell=True)
    p4 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam04), shell=True)
    p5 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam05), shell=True)
    p6 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam06), shell=True)
    p7 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam07), shell=True)
    p8 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam08), shell=True)
    p9 = subprocess.Popen("samtools view -H {0} | cat - {1}.tmp | samtools view -bS - | samtools sort -f - {1}".format(bam, bam09), shell=True)

    p1.wait()
    p2.wait()
    p3.wait()
    p4.wait()
    p5.wait()
    p6.wait()
    p7.wait()
    p8.wait()
    p9.wait()


    p1 = subprocess.Popen("rm {}.tmp".format(bam01), shell=True)
    p2 = subprocess.Popen("rm {}.tmp".format(bam02), shell=True)
    p3 = subprocess.Popen("rm {}.tmp".format(bam03), shell=True)
    p4 = subprocess.Popen("rm {}.tmp".format(bam04), shell=True)
    p5 = subprocess.Popen("rm {}.tmp".format(bam05), shell=True)
    p6 = subprocess.Popen("rm {}.tmp".format(bam06), shell=True)
    p7 = subprocess.Popen("rm {}.tmp".format(bam07), shell=True)
    p8 = subprocess.Popen("rm {}.tmp".format(bam08), shell=True)
    p9 = subprocess.Popen("rm {}.tmp".format(bam09), shell=True)
    p10 = subprocess.Popen("rm {0}.shuffled.sam".format(bam), shell=True)


    p1.wait()
    p2.wait()
    p3.wait()
    p4.wait()
    p5.wait()
    p6.wait()
    p7.wait()
    p8.wait()
    p9.wait()
    p10.wait()

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




    args = parser.parse_args()

    bam01, bam02 = pre_process_bam(args.bam, args.bam01, args.bam02, args.bam03,
                                   args.bam04, args.bam05, args.bam06, args.bam07,
                                   args.bam08, args.bam09)
