####
#
#	Count tags to regions script. Not optimized for parallel
#	environments. Produces a genic counts file. 
#
#	ppliu
#
#####

import csv
from collections import defaultdict, namedtuple
import itertools
from optparse import OptionParser
from subprocess import Popen, PIPE
import os


from clipper.src.readsToWiggle import readsToWiggle_pysam
import numpy as np
import pysam


				
def union(*dicts):
	
	""" unions a list of dicts into one dict """
	return dict(itertools.chain(*map(lambda dct: list(dct.items()), dicts)))
   
def count_to_regions(basedir, species):

	"""

    Gets all genic regions, returns a two dicts, a regions dict and a genes dict
	for all the genes

	returns genes
	
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

	genes = defaultdict(lambda : {'regions' : [],
				      'start' : np.inf,
				      'stop'  : np.NINF,
				      'chr'   : None,
				      'strand': None,
				      'frea'  : None,
				      'raw_count' : None,
				      'gene_id' : None}
				      )
	
	for chr in chrs:
		regions_file = os.path.join(basedir, "genic_regions_"+species+".chr"+chr)
		with open(regions_file, 'r') as gene_file:
			for line in csv.reader(gene_file, delimiter="\t"):
				
				chromosome, start, stop, strand, ensembl_id, frea_annot = line

				strand = "+" if int(strand) == 0 else "-"
				
				genes[ensembl_id]['regions'].append((start, stop))
	
				#get the minimal start, and maximal stop for each gene
				genes[ensembl_id]['start'] = min(int(start), genes[ensembl_id]["start"])
				genes[ensembl_id]['stop'] = max(int(stop), genes[ensembl_id]["stop"])
				genes[ensembl_id]["chr"] = chromosome
				genes[ensembl_id]["strand"] = strand
				genes[ensembl_id]["frea"] = frea_annot
				genes[ensembl_id]["raw_count"] = 0
				genes[ensembl_id]['gene_id'] = ensembl_id
	
	return genes

def count_gene(bam_file, gene, flip):
	
	"""
	
	get read counts for genic regions in the gene specified by annotation in passed value 'keys'

	"""
	count = namedtuple('count', ['gene_count', 'region_count'])
    
	region_counts = {}

	try:
		# fetch reads from bam file for the gene referenced by keys (Ensembl ID)
		subset_reads = bam_file.fetch(reference = gene['chr'],
					      start = gene["start"],
					      end = gene["stop"])
		
	except:
		raise Exception("could not fetch reads. check if bam is indexed %s:%s-%s" % (gene['chr'], gene['start'], gene['stop']))

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
	for region_start, region_stop in gene['regions']:
		
		start = int(region_start) - gene["start"]
		stop  = int(region_stop)  - gene["start"]
		
		gene_sum += sum(wig[start:stop])
		
		region_counts[gene['gene_id'] + str(start + gene["start"]) + str(stop + gene["start"])] = sum(wig[start:stop])
	
	return {region : count(gene_sum, region_sum) for region, region_sum in region_counts.items()}

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
	basedir = os.path.join(basedir, "ppliu/genic_regions/"+species+"/")
	genes = count_to_regions(basedir, species)
	
	region_counts = []
	for gene in genes.values():
		region_counts.append(count_gene(bam_file, gene, flip)) 
	
	region_counts = union(region_counts)
	with open(out_file, 'w') as out_file:
		for line in csv.reader(genic_order_file, delimiter="\t"):
			chrom, start, stop, strand, ensembl_id, frea_annot = line
			
			if ensembl_id+start+stop in region_counts:
				out_file.write("\t".join([

										str(chrom), str(start), str(stop), 
										str(region_counts[ensembl_id+start+stop]), 
										str(gene_counts[ensembl_id]), 
										str(strand), str(ensembl_id), str(frea_annot), "\n"
										])
							)

if __name__ == "__main__":
	# detect between oolite and triton hosts
	host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()

	if "optiputer" in host or "compute" in host:
		basedir = "/nas/nas0/"
	elif "tcc" in host or "triton" in host:
		basedir = "/projects/"
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
