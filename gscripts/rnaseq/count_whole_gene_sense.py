####
#
#	Q&D script to count reads to whole gene space, and count to sense or
#	antisense.
#
#	ppliu
#
####

# import dependencies
import pysam
import pybedtools
from optparse import OptionParser
# initialize global variables
gene_info = dict()			# holds info about gene itself

# main function
def main():
	
	# gather command line option values
	parser = OptionParser()
	parser.add_option("-b", "--bam_file", dest="bam_path")
	parser.add_option("-s", "--species", dest="species")

	# assign option values to variables
	(options, args) = parser.parse_args()
	bam_path = options.bam_path
	species = options.species
		
	# open alignment file
	bam_file = pysam.Samfile(bam_path, "rb")



	get_genic_regions(species)


	antisense_total=0
	sense_total=0
	total_reads=0

	for keys in gene_info:

		if (keys not in gene_info):
			print keys
			print "not here"

		chr_index = gene_info[keys]["chr"]
		start = gene_info[keys]["start"]
		stop = gene_info[keys]["stop"]
		strand = gene_info[keys]["strand"]

		# fetch reads from region in alignment file
		iter = bam_file.fetch( chr_index, int(start), int(stop) )

		# track count of sense and antisense reads
		antisense_count=0
		sense_count=0

		for read in iter: 

			# if sense strand
			if ( (read.is_reverse and strand == "-") or (not read.is_reverse and strand =="+")):
				
				sense_count+=1

			elif ( (read.is_reverse and strand == "+") or (not read.is_reverse and strand =="-")):

				antisense_count+=1

				
			total_reads+=1

		antisense_total+=antisense_count
		sense_total+=sense_count

		if (total_reads > 0):
			antisense_percentage = float(antisense_total)/float(total_reads)
			sense_percentage = float(sense_total)/float(total_reads)

		else:
			antisense_percentage = 0
			sense_percentage = 0

		print str(keys)+"\tA:"+str(antisense_total)+"\tS:"+str(sense_total)+"\tT:"+str(total_reads)+"\tA:"+str(antisense_percentage)+"\tS:"+str(sense_percentage)

	print "sense count: "+str(sense_total)
	print "anti-sense count: "+str(antisense_total)
	print "sense percentage: "+str(sense_percentage)
	print "anti-sense percentage: "+str(antisense_percentage)

	return

# function to fill dictionary with annotated genic regions
def get_genic_regions(species):

	# for each species, define the chromosome indexes
	if species == "hg19":
		chrs = map(str,range(1,23)) #1-22
		chrs.append("X")
		chrs.append("Y")

	elif species == "mm9":
		chrs = map(str, range(1,20))
		chrs.append("X")
		chrs.append("Y")

	elif species == "ce6":
		chrs = ("I", "II", "III", "IV", "V", "X")

	# open the genic regions annotation file
	for chr in chrs:
		regions_file = "/projects/ppliu/genic_regions/"+species+"/genic_regions_"+species+".chr"+chr
		f = open(regions_file, 'r')
		for line in f:

			line = line.strip()
			chromosome, start, stop, strand, ensembl_id, frea_annot = line.strip().split("\t")
	
			if ( int(strand) == 1):
				strand = '-'
			elif( int(strand) == 0):
				strand = '+'

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

if __name__ == '__main__': 
	
	main()
