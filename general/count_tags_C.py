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
import pickle
import os
from optparse import OptionParser
from subprocess import Popen, PIPE
import sys
import random
import re
import time

from bx.bbi.bigwig_file import BigWigFile
from clipper.src.peaks import readsToWiggle_pysam
from numpy import *
import pybedtools
import pysam

# get read counts for genic regions in the gene specified by annotation in passed value 'keys'
def get_wig(key, bam_file, gene_info, regions_info, flip):

	if key not in gene_info:
		print keys
		print "not here"

	try:
		# fetch reads from bam file for the gene referenced by keys (Ensembl ID)
		subset_reads = bam_file.fetch(reference = gene_info[keys]['chr'],
					      start = int(gene_info[keys]["start"]),
					      end = int(gene_info[keys]["stop"]))
		
	except:
		raise Exception("could not fetch reads. check if bam is indexed."+gene_info[key]['chr']+" "+ gene_info[key]['start']+" "+ gene_info[key]['stop'])

	try:
		# determine strand to keep based on flip option
		keep_strand = gene_info[key]["strand"]
		if str(flip) == "flip":
			if str(keep_strand) == '-':
				keep_strand = '+'
			elif str(keep_strand) == '+':
				keep_strand = '-'

		elif str(flip) == "both":
			keep_strand = 0;

		wig, jxns, nr_counts, read_lengths, reads = readsToWiggle_pysam(subset_reads,
										int(gene_info[key]["start"]),
										int(gene_info[key]["stop"]),
										keep_strand,
										'center',
										True)

	except Exception as e:
		print e

	region_sum = 0
	while len(regions_info[key]) >= 2:
		try:
			start = int(regions_info[key].pop(0)) - int(gene_info[key]["start"] )
			stop = int(regions_info[key].pop(0)) - int(gene_info[key]["start"] )
		except:
			print "no regions information"

		sum = 0
		for x in range(start, stop):
			sum+=wig[x]	
		
		region_sum+=sum
		gene_info[key+str(int(start)+int(gene_info[keys]["start"] ))+str(int(stop)+int(gene_info[key]["start"]) )] =  str(sum)

	gene_info[key]["region_sum"] = str(region_sum)
	

def count_to_regions(basedir, species):

	"""

	Gets all genic regions, returns a two dicts, a regions dict and a genes dict
	for all the genes

	returns gene_info, regions_info
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

	gene_info = defaultdict(list)
	regions_info = defaultdict(list)
	
	for chr in chrs:
		regions_file = basedir+"/ppliu/genic_regions/"+species+"/genic_regions_"+species+".chr"+chr
		with open(regions_file, 'r') as gene_file:
			for line in csv.reader(gene_file, delimiter="\t"):

				chromosome, start, stop, strand, ensembl_id, frea_annot = line
	
				strand = "+" if int(strand) == 0 else "-"
			
				regions_info[ensembl_id].append(start)
				regions_info[ensembl_id].append(stop)

				if (ensembl_id in gene_info):
					#get the minimal start, and maximal stop for each gene
					gene_info[ensembl_id]['start'] = min(int(start),
									       int(gene_info[ensembl_id]["start"]))
						
					gene_info[ensembl_id]['stop'] = max(int(stop),
									     int(gene_info[ensembl_id]["stop"]))
				else:
					gene_info[ensembl_id]={}
					gene_info[ensembl_id]["chr"] = chromosome
					gene_info[ensembl_id]["start"] = start
					gene_info[ensembl_id]["stop"] = stop
					gene_info[ensembl_id]["strand"] = strand
					gene_info[ensembl_id]["frea"] = frea_annot
					gene_info[ensembl_id]["raw_count"] = 0
	
	return gene_info, regions_info

def main(genic_order_file):
	
	region_lines = [get_wig(key) for key in gene_info.keys()]

	for lines in genic_order_file:

		lines = lines.strip()
		chr, start, stop, strand, ensembl_id, frea_annot = lines.split("\t")
		
		if ensembl_id+start+stop in gene_info:
			out_file.write( chr+"\t"+start+"\t"+stop+"\t"+gene_info[ensembl_id+start+stop]+"\t"+gene_info[ensembl_id]["region_sum"]+"\t"+strand+"\t"+ensembl_id+"\t"+frea_annot+"\n")

if __name__ == "__main__":

	# detect between oolite and triton hosts
	host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
	if "optiputer" in host or "compute" in host:
		basedir = "/nas/nas0"
	elif "tcc" in host or "triton" in host:
		basedir = "/projects"
	else:
		print "Not in the correct location"
		raise Exception

	# gather command line option values
	parser = OptionParser()
	parser.add_option("-s", "--species", dest="species")
	parser.add_option("-b", "--bam_file", dest="bam_path")
	parser.add_option("-f", "--flip", dest="flip")
	parser.add_option("-o", "--out_file", dest="out_file")
	
	# assign parameters to variables
	(options,args) = parser.parse_args()
	flip = options.flip
	
	# open output file for writing 
	out_file = open(options.out_file, 'a')
	
	# open properly ordered genic order file for reading
	genic_order_file = open(basedir+"/ppliu/genic_counts_orders/"+species+".order", 'r')
	
	# open bam file for reading
	bam_file = pysam.Samfile(options.bam_path, 'rb')
	
	# create dictionary data structures 
	gene_info, regions_info = count_to_regions(basedir, options.species)
	main()
	
