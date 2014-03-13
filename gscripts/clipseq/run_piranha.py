#Given a bam file compares clipper and piranah to same outputs
#will want to better paramaterize this later...
from optparse import OptionParser
import os
import subprocess
import sys

import pybedtools

parser = OptionParser()
parser.add_option("-b", "--bam", dest="bam", help="bam file to run peak finding algorithm on")
parser.add_option("-p", "--p_value", dest="p_value", default=.05, help="p_value cutoff for peaks")
parser.add_option("-s", "--bin_size", dest="bin_size", default=200)
parser.add_option("-o", "--out_file", dest="out_file")



# assign parameters to variables
(options,args) = parser.parse_args()
                                                               
#format reads for piranah, this is a bitch
s = subprocess.Popen("bamToBed -i %s | sort -k1,1 -k3n,3 -k2,2n > %s_tmp.bed" % (options.bam, options.bam), shell=True)
s.wait()
s = subprocess.Popen("cut -f 1,2,3 %s_tmp.bed > %s_input.bed" % (options.bam, options.bam), shell=True) 
s.wait()
s = subprocess.Popen("Piranha -b %s %s_input.bed > %s_results.bed" % (options.bin_size, options.bam, options.bam), shell=True) 
s.wait()
s = subprocess.Popen("cat %s_results.bed | awk '{if($7 < %s ) print $0}' | awk 'BEGIN {OFS=\"\t\"} {print $1, $2, $3, $4, $7, $6}' > %s" % (options.bam, options.p_value, options.out_file), shell=True)
s.wait()
#cleanup

os.remove("%s_tmp.bed" % (options.bam))
os.remove("%s_input.bed" % (options.bam))
os.remove("%s_results.bed" % (options.bam))
