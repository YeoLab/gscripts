import sys
from conserved_structure import RNApair
import pybedtools

fi = "/nas3/lovci/projects/structure/hg19/all.RNAlinks.sorted.bed"
out1 = open("out1", 'w')
out2 = open("out2", 'w')
for line in pybedtools.BedTool(fi):
    pair = RNApair(line)
    out1.write("\t".join([pair.chromosome, str(pair.qStart), str(pair.qStop), pair.name]) + "\n")
    out2.write("\t".join([pair.chromosome, str(pair.tStart), str(pair.tStop), pair.name]) + "\n") 
