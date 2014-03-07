import pybedtools
import argparse

OFFSET = 14

def five_prime(interval):

    #There could be an off by one bug here, need to think about it a bit
    if interval.strand == "+":
        interval.start = interval.start + OFFSET
        interval.end = interval.start + 1
    else:
        interval.end = interval.end - OFFSET 
        interval.start = interval.end - 1

    return interval

parser = argparse.ArgumentParser(description="Adjusts bam files (or bed files) to correct position for riboseq, outputs bed file of corrected positions")
parser.add_argument("--bam", help="bam file to adjust", required=True)
parser.add_argument("--out", help="output file (bed format)", required=True)
args = parser.parse_args()

pybedtools.BedTool(args.bam).bam_to_bed().each(five_prime).sort().saveas(args.out) 
