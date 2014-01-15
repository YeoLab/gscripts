import pybedtools
import argparse

def adjust_start(feature):
    feature.start = feature.start + 13
    feature.end = feature.end + 13
    return feature

def five_prime(feature):
    if feature.strand == "+":
        feature.end = feature.start
    else:
        feature.start = feature.end
    return feature

parser = argparse.ArgumentParser(description="Adjusts bam files (or bed files) to correct position for riboseq, outputs bed file of corrected positions")
parser.add_argument("--bam", help="bam file to adjust", required=True)
parser.add_argument("--out", help="output file (bed format)", required=True)
args = parser.parse_args()

pybedtools.BedTool(args.bam).bam_to_bed().each(five_prime).each(adjust_start).sort().saveas(args.out) 
