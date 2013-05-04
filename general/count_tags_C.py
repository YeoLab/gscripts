####
#
#	Count tags to regions script. Not optimized for parallel
#	environments. Produces a genic counts file. 
#
#	ppliu
#
#####

import csv
from collections import defaultdict
from optparse import OptionParser
from subprocess import Popen, PIPE

from clipper.src.call_peak import readsToWiggle_pysam
import numpy as np
import pysam

def count_reads(bam_file, genes, regions, flip):
	
	"""
	
	get read counts for genic regions in the gene specified by annotation in passed value 'keys'

	"""
	
	gene_counts = {}
	region_counts = {}
	for key, gene in genes.items():

		try:
			# fetch reads from bam file for the gene referenced by keys (Ensembl ID)
			subset_reads = bam_file.fetch(reference = gene['chr'],
						      start = gene["start"],
						      end = gene["stop"])
			
		except:
			raise Exception("could not fetch reads. check if bam is indexed." + 
						gene['chr']+ " " + gene['start']+" "+ gene['stop'])
	
		# determine strand to keep based on flip option
		keep_strand = gene["strand"]
		if str(flip) == "flip":
			if str(keep_strand) == '-':
				keep_strand = '+'
			elif str(keep_strand) == '+':
				keep_strand = '-'

		elif str(flip) == "both":
			keep_strand = 0;

		wig, jxns, nr_counts, read_lengths, reads = readsToWiggle_pysam(
																	subset_reads,
																	int(gene["start"]),
																	int(gene["stop"]),
																	keep_strand,
																	'center',
																	True)


	
		gene_sum = 0
		for region_start, region_stop in regions[key]:
			
			start = int(region_start) - gene["start"]
			stop  = int(region_stop)  - gene["start"]
			
			gene_sum += sum(wig[start:stop])
			
			region_counts[key + str(start + gene["start"]) + str(stop + gene["start"])] =  sum(wig[start:stop])
	
		gene_counts[key] = gene_sum
	return gene_counts, region_counts

def count_to_regions(basedir, species):

	"""

	Gets all genic regions, returns a two dicts, a regions dict and a genes dict
	for all the genes

	returns genes, regions
	
	"""

	if species == "hg19" or species == "hg18":
		chrs = map(str,range(1,23)) #1-22
		chrs.append("X")
		chrs.append("Y")

	elif species == "mm9":
		chrs = map(str, range(1,20))
		chrs.append("X")
		chrs.append("Y")
	elif species == "ce6":
		chrs = ("I", "II", "III", "IV", "V", "X")

	genes = defaultdict(list)
	regions = defaultdict(list)
	
	for chr in chrs:
		regions_file = basedir+"/ppliu/genic_regions/"+species+"/genic_regions_"+species+".chr"+chr
		with open(regions_file, 'r') as gene_file:
			for line in csv.reader(gene_file, delimiter="\t"):

				chromosome, start, stop, strand, ensembl_id, frea_annot = line
	
				strand = "+" if int(strand) == 0 else "-"
			
				regions[ensembl_id].append((start, stop))

				if (ensembl_id in genes):
					#get the minimal start, and maximal stop for each gene
					genes[ensembl_id]['start'] = min(int(start),
									       int(genes[ensembl_id]["start"]))
						
					genes[ensembl_id]['stop'] = max(int(stop),
									     int(genes[ensembl_id]["stop"]))
				else:
					genes[ensembl_id]={}
					genes[ensembl_id]["chr"] = chromosome
					genes[ensembl_id]["start"] = int(start)
					genes[ensembl_id]["stop"] = int(stop)
					genes[ensembl_id]["strand"] = strand
					genes[ensembl_id]["frea"] = frea_annot
					genes[ensembl_id]["raw_count"] = 0
	
	return genes, regions

def count_tags(basedir, species, bam_file, flip, out_file):
	
	"""
		Main function counts tags and ouptouts counts to outfile
		
		basedir - root where files are stored (generally /nas/nas0
		species - species to count for, hg19, hg18, mm9
		flip - if true flip the reads
		out_file - output file 
		
	"""
	# open properly ordered genic order file for reading
	genic_order_file = open(basedir+"/ppliu/genic_counts_orders/"+species+".order", 'r')
	
	# open bam file for reading
	bam_file = pysam.Samfile(bam_file, 'rb')
	
	# create dictionary data structures 
	genes, regions = count_to_regions(basedir, species)
	
	gene_counts, region_counts = count_reads(bam_file, genes, regions, flip) 
	
	with open(out_file, 'w') as out_file:
		for line in csv.reader(genic_order_file, delimiter="\t"):
			chrom, start, stop, strand, ensembl_id, frea_annot = line
			
			if ensembl_id+start+stop in region_counts:
				out_file.write("\t".join([
										chrom, start, stop, 
										region_counts[ensembl_id+start+stop], 
										gene_counts[ensembl_id], 
										strand, ensembl_id, frea_annot, "\n"
										])
							)

if __name__ == "__main__":
	# detect between oolite and triton hosts
	host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()

	if "optiputer" in host or "compute" in host:
		basedir = "/nas/nas0"
	elif "tcc" in host or "triton" in host:
		basedir = "/projects"
	else:
		raise Exception("Not in the correct location, current host: %s" % (host))

	# gather command line option values
	parser = OptionParser()
	parser.add_option("-s", "--species", dest="species")
	parser.add_option("-b", "--bam_file", dest="bam_path")
	parser.add_option("-f", "--flip", dest="flip", action="store_true", help="Flip reads", default=False)
	parser.add_option("-o", "--out_file", dest="out_file")
	
	# assign parameters to variables
	(options,args) = parser.parse_args()
	count_tags(basedir, options.species, options.bam_path, options.flip, options.out_file)
