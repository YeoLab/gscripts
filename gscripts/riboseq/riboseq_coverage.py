import pybedtools
import argparse

OFFSET = 14

def five_prime(feature):
    if feature.strand == "+":
        feature.start = feature.start + OFFSET
        feature.end = feature.start
    else:
        feature.end = feature.end - OFFSET
        feature.start = feature.end
    return feature

parser = argparse.ArgumentParser(description="Adjusts bam files (or bed files) to correct position for riboseq, outputs bed file of corrected positions")
parser.add_argument("--bam", help="bam file to adjust", required=True)
parser.add_argument("--out", help="output file (bed format)", required=True)
args = parser.parse_args()

pybedtools.BedTool(args.bam).bam_to_bed().each(five_prime).sort().saveas(args.out) 
