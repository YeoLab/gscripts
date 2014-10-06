#Given a bam file compares clipper and piranah to same outputs
#will want to better paramaterize this later...
import subprocess
import sys
import pybedtools
from optparse import OptionParser
import os
import shutil

parser = OptionParser()
parser.add_option("-b", "--bam", dest="bam", help="bam file to run peak finding algorithm on")
parser.add_option("-s", "--species", dest="species")
parser.add_option("-o", "--out_file", help="out file")
# assign parameters to variables
(options,args) = parser.parse_args()

s = subprocess.call("bamToBed -i %s > %s_tmp.bed" % (options.bam, options.bam), shell=True)
if options.species == "hg19": 
    s = subprocess.call("perl ~/gscripts/gscripts/clipseq/kasey_peak_calling/allocate_reads_to_gene_hg19_fromBED.pl hg19 %s_tmp.bed" % (options.bam), shell = True) 
    
elif options.species == "mm9":
    s = subprocess.call("perl ~/gscripts/gscripts/clipseq/kasey_peak_calling/allocate_reads_to_gene_mm9_fromBED.pl mm9 %s_tmp.bed" % (options.bam), shell = True) 
else:
    raise "species not supported"



subprocess.call("perl ~/gscripts/gscripts/clipseq/kasey_peak_calling/analyze_clip_data_complete %s %s_tmp.bed.ingenes.BED -force -shift 0 -trim 0 -window_min 50 -window_max 150 -gap 20 -mRNA 0 -pval 0.0001 0" % (options.species, options.bam), shell = True)

shutil.copy("%s_tmp.bed_data/%s_tmp.bed_notrim_s0_50-150_g20_corrected_premRNA_withrepeats_clusters/%s_tmp.bed_notrim_ingenes_clusters_%s50.bed" % (options.bam, options.bam, options.bam, options.species), options.out_file) 
