####
#
#	Count tags to regions script. Not optimized for parallel
#	environments. Produces a genic counts file. 
#
#	ppliu
#
#####

print "import"

# import dependencies
import pybedtools
import pickle
import random
from numpy import *
import subprocess
from bx.bbi.bigwig_file import BigWigFile
from subprocess import Popen, PIPE
import re
import time
import pysam
#from seqTools import *
import sys
import os
from optparse import OptionParser
from clipper.src.peaks import readsToWiggle_pysam

print "done import"

print "detect host"
# detect between oolite and triton hosts
host = Popen(["hostname"], stdout=PIPE).communicate()[0].strip()
if "optiputer" in host or "compute" in host:
	basedir = "/nas/nas0"
elif "tcc" in host or "triton" in host:
	basedir = "/projects"
else:
	print "Where am I?"
	raise Exception

print "done detect host"

print "parse options"
# gather command line option values
parser = OptionParser()
parser.add_option("-s", "--species", dest="species")
parser.add_option("-b", "--bam_file", dest="bam_path")
parser.add_option("-f", "--flip", dest="flip")
parser.add_option("-o", "--out_file", dest="out_file")
print "done parse options"


 # assign parameters to variables
(options,args) = parser.parse_args()
species = options.species
bam_path = options.bam_path
flip = options.flip

# open output file for writing 
out_file = open(options.out_file, 'a')

# open properly ordered genic order file for reading
genic_order_file = open(basedir+"/ppliu/genic_counts_orders/"+species+".order", 'r')

print "open bam"
# open bam file for reading
bam_file = pysam.Samfile(bam_path, 'rb')
print "done open bam"

# create dictionary data structures 
gene_info = dict()
regions_info = dict()


# get read counts for genic regions in the gene specified by annotation in passed value 'keys'
def get_wig(keys):

	if (keys not in gene_info):
		print keys
		print "not here"

	try:

		# fetch reads from bam file for the gene referenced by keys (Ensembl ID)
		subset_reads = bam_file.fetch(reference = gene_info[keys]['chr'], start = int(gene_info[keys]["start"]), end = int(gene_info[keys]["stop"]))
		
	except:
		print "could not fetch reads. check if bam is indexed."+gene_info[keys]['chr']+" "+ gene_info[keys]['start']+" "+ gene_info[keys]['stop']

	try:
		# determine strand to keep based on flip option
		keep_strand = gene_info[keys]["strand"]
		if (str(flip) == "flip"):
			if (str(keep_strand) == '-'):
				keep_strand = '+'
			elif (str(keep_strand) == '+'):
				keep_strand = '-'

		elif (str(flip) == "both"):
			keep_strand = 0;

		wig, jxns, nr_counts, read_lengths, reads = readsToWiggle_pysam(subset_reads, int(gene_info[keys]["start"]), int(gene_info[keys]["stop"]), keep_strand, 'center', True)

	except Exception as e:
		print e

	region_sum=0
	while (len(regions_info[keys]) >=2 ):

		try:
			start = int(regions_info[keys].pop(0)) - int(gene_info[keys]["start"] )
			stop = int(regions_info[keys].pop(0)) - int(gene_info[keys]["start"] )
		except:
			print "no regions information"

		sum = 0
		for x in range(start, stop):
			sum+=wig[x]	
		
		region_sum+=sum
		gene_info[keys+str(int(start)+int(gene_info[keys]["start"] ))+str(int(stop)+int(gene_info[keys]["start"]) )] =  str(sum)

	gene_info[keys]["region_sum"] = str(region_sum)
	

def count_to_regions(species):

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
	for chr in chrs:
		regions_file = basedir+"/ppliu/genic_regions/"+species+"/genic_regions_"+species+".chr"+chr
		f = open(regions_file, 'r')
		for line in f:

			line = line.strip()
			chromosome, start, stop, strand, ensembl_id, frea_annot = line.strip().split("\t")
	

			if ( int(strand) == 1):
				strand = '-'
			elif( int(strand) == 0):
				strand = '+'
		
			if (ensembl_id not in gene_info):
				regions_info[ensembl_id] = []

			regions_info[ensembl_id].append(start)
			regions_info[ensembl_id].append(stop)



			if (ensembl_id in gene_info):
				if (int(start) < int(gene_info[ensembl_id]["start"])):
					gene_info[ensembl_id]["start"] = start
				if (int(stop) > int(gene_info[ensembl_id]["stop"])):
					gene_info[ensembl_id]["stop"] = stop

			else:
				gene_info[ensembl_id]={}
				gene_info[ensembl_id]["chr"] = chromosome
				gene_info[ensembl_id]["start"] = start
				gene_info[ensembl_id]["stop"] = stop
				gene_info[ensembl_id]["strand"] = strand
				gene_info[ensembl_id]["frea"] = frea_annot
				gene_info[ensembl_id]["raw_count"] = 0
		f.close()

def main():
	
	region_lines = [get_wig(key) for key in gene_info.keys()]

	for lines in genic_order_file:

		lines = lines.strip()
		chr, start, stop, strand, ensembl_id, frea_annot = lines.split("\t")
		
		if ensembl_id+start+stop in gene_info:
			out_file.write( chr+"\t"+start+"\t"+stop+"\t"+gene_info[ensembl_id+start+stop]+"\t"+gene_info[ensembl_id]["region_sum"]+"\t"+strand+"\t"+ensembl_id+"\t"+frea_annot+"\n")

if __name__ == "__main__":

	count_to_regions(species)
	main()
	
