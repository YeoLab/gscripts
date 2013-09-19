import pybedtools
import argparse

parser = argparse.ArgumentParser(description="Calculates Non-redundant fraction of reads in a bam file")
parser.add_argument("--bam", help="bam file to calculate NRF from", required=True)
parser.add_argument("--genome", help="Genome sizes", required=True)
args = parser.parse_args()

coverage = pybedtools.BedTool(args.bam)
NRF = coverage.genome_coverage(g=args.genome, **{'5':True})
reads_mapped = 0.0
total_locations_mapped = 0.0
locations_with_one = 0.0
#makes historgram
#chr number count
for line in str(NRF).strip().split("\n"):
        line = line.strip().split()
        number_of_reads = int(line[1])
	locations_mapped = int(line[2])
	if number_of_reads != 0:
		total_locations_mapped += locations_mapped
		reads_mapped += number_of_reads * locations_mapped
		if number_of_reads == 1:
			locations_with_one += locations_mapped
print "NRF", "PCB"
print (total_locations_mapped / reads_mapped), (locations_with_one / total_locations_mapped)
