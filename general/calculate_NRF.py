import pybedtools
import argparse

parser = argparse.ArgumentParser(description="Calculates Non-redundant fraction of reads in a bam file")
parser.add_argument("--bam", help="bam file to calculate NRF from", required=True)
parser.add_argument("--genome", help="Genome sizes", required=True)
args = parser.parse_args()

coverage = pybedtools.BedTool(args.bam)
NRF = coverage.genome_coverage(g=args.genome, **{'5':True})
reads_mapped = 0.0
locations_mapped = 0.0
for line in str(NRF).strip().split("\n"):
        line = line.strip().split()
        if int(line[1]) != 0:
            locations_mapped += int(line[2])
            reads_mapped += int(line[1]) * int(line[2])

print "NRF"
print locations_mapped / reads_mapped
