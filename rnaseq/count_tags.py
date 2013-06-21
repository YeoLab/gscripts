####
#
#	Count tags to regions script. 
#	environments. Produces a genic counts file. 
#
#	ppliu,G gabriel Pratt
#
#####

import csv
from collections import defaultdict, namedtuple
import itertools
from optparse import OptionParser
from subprocess import Popen, PIPE
import multiprocessing 
import os


from clipper.src.readsToWiggle import readsToWiggle_pysam
import numpy as np
import pysam

count = namedtuple('count', ['gene_count', 'region_count'])
				
def union(*dicts):
	
	""" unions a list of dicts into one dict """
	return dict(itertools.chain(*map(lambda dct: list(dct.items()), dicts)))
   
def count_to_regions(annotation):

	"""

    Gets all genic regions, returns a two dicts, a regions dict and a genes dict
	for all the genes

	returns genes
	
	"""

	genes = defaultdict(lambda : {'regions' : [],
				      'start' : np.inf,
				      'stop'  : np.NINF,
				      'chr'   : None,
				      'strand': None,
				      'frea'  : None,
				      'raw_count' : None,
				      'gene_id' : None}
				      )
	

	with open(annotation, 'r') as gene_file:
		for line in csv.reader(gene_file, delimiter="\t"):
			
			chromosome, start, stop, name, score, strand, exon_number = line
			
			strand = "+" if int(strand) == 0 else "-"
			
			genes[ensembl_id]['regions'].append((start, stop))
			
			#get the minimal start, and maximal stop for each gene
			genes[ensembl_id]['start'] = min(int(start), genes[ensembl_id]["start"])
			genes[ensembl_id]['stop'] = max(int(stop), genes[ensembl_id]["stop"])
			genes[ensembl_id]["chr"] = chromosome
			genes[ensembl_id]["strand"] = strand
			genes[ensembl_id]["frea"] = exon_number
			genes[ensembl_id]["raw_count"] = 0
			genes[ensembl_id]['gene_id'] = ensembl_id
	
	return genes

def count_gene(bam_file, gene, flip):
	
	"""
	
	get read counts for genic regions in the gene specified by annotation in passed value 'keys'

	bam_file - pysam bam file
	
	"""
	
	
    
	region_counts = {}


	#try:
	bam_file = pysam.Samfile(bam_file, 'rb')
	# fetch reads from bam file for the gene referenced by keys (Ensembl ID)
	subset_reads = bam_file.fetch(reference = gene['chr'],
				      start = int(gene["start"]),
				      end = int(gene["stop"]))
		
	#except:
	#	raise Exception("could not fetch reads. check if bam is indexed %s:%s-%s" % (gene['chr'], gene['start'], gene['stop']))


	# determine strand to keep based on flip option
	keep_strand = gene["strand"]
	if str(flip) == "flip":
		if str(keep_strand) == '-':
			keep_strand = '+'
		elif str(keep_strand) == '+':
			keep_strand = '-'

	elif str(flip) == "both":
		keep_strand = 0;

	wig, jxns, nr_counts, read_lengths, reads = readsToWiggle_pysam(subset_reads,
									int(gene["start"]),
									int(gene["stop"]),
									keep_strand,
									'center', True)

	gene_sum = 0
	for region_start, region_stop in gene['regions']:
		
		start = int(region_start) - gene["start"]
		stop  = int(region_stop)  - gene["start"]
		
		gene_sum += sum(wig[start:stop])
		
		region_counts[gene['gene_id'] + str(start + gene["start"]) + str(stop + gene["start"])] = sum(wig[start:stop])

	bam_file.close()
	return [(region, count(gene_sum, region_sum)) for region, region_sum in region_counts.items()]

def func_star(varables):
	""" covert f([1,2]) to f(1,2) """
	return count_gene(*varables)

def count_tags(bam_file, filp, out_file, annotation, num_cpu = "autodetect", ):
	
	"""
		Main function counts tags and ouptouts counts to outfile
		
		bam_file - bam file to process
		flip - if true flip the reads
		out_file - output file 
		num_cpu - number of cpus to use in parallel processing, autodetect uses all avaiable cpus
		annotation - path to annotation file, format is bed 6 + 1 where 1 is exon information (may swap this out later)
	"""
	
	if num_cpu == 'autodetect':
		num_cpu = multiprocessing.cpu_count()
		
	
	if not os.path.exists(bam_file):
		raise Exception("bam file %s does not exist" % (bam_file))
				
	genes = count_to_regions(annotation)
	
	#region_counts =  [count_gene(bam_file, gene, flip) for gene in genes.values()]

	pool = multiprocessing.Pool(int(num_cpu))
			
	region_counts = pool.map(func_star, [(bam_file, gene, flip) for gene in genes.values()], chunksize=50)	
	#region_counts = [job.get(timeout=360) for job in jobs]
	
	region_counts = dict(itertools.chain(*region_counts))
	
	with open(out_file, 'w') as out_file:
		for chrom, start, stop, strand, ensembl_id, frea_annot in csv.reader(genic_order_file, delimiter="\t"):
			if ensembl_id+start+stop in region_counts:
				out_file.write("\t".join([str(chrom), str(start), str(stop), 
							  str(region_counts[ensembl_id+start+stop].region_count), 
							  str(region_counts[ensembl_id+start+stop].gene_count), 
							  str(strand), str(ensembl_id), str(frea_annot), "\n"
							  ]))

if __name__ == "__main__":

	# gather command line option values
	parser = OptionParser()
	parser.add_option("-b", "--bam_file", dest="bam_path")
	parser.add_option("-f", "--flip", dest="flip", help="Flip reads to flip type : flip, for non-strand specific: both", default=False)
	parser.add_option("-o", "--out_file", dest="out_file")
	parser.add_option("--processors", dest="np", default="autodetect",
			  help="Number of processors to use. Default: All processors on machine",
			  type="str", metavar="NP")
	parser.add_option("--annotation_file", dest="annotation", help="annotation to count tags from, generated from gtfutils")
	
	# assign parameters to variables
	(options,args) = parser.parse_args()
	count_tags(options.bam_path, options.flip,
		   options.out_file, num_cpu = options.np,
		   annotation = options.annotation)
