#This script requieres at least a BED-4 (a unique name) to run properly
#need to figure out how to handle arbitray bed files...
#This is actually to slow to run efficently 
import sys
import pybedtools


def convert_to_mRNA(feature):
    pass

def makeStarts(feature):
    feature.end = feature.start + 1
    return feature

def makeEnds(feature):
    feature.start = feature.end - 1
    return feature

def convert_to_mRNA(feature):
    feature.chrom = feature[10]
        
    if feature.strand == "+":
        #start of exon location in mRNA + distance from exon start in genomic rna
        start = int(feature[11]) + (int(feature.start) - int(feature[7]))
        end   = int(feature[11]) + (int(feature.end) - int(feature[7]))
        
    else: #strand = -
        
        start = int(feature[11]) + (int(feature[8]) - int(feature.end))
        end   = int(feature[11]) + (int(feature[8]) - int(feature.start))
        
    feature.start = start
    feature.end   = end
    return feature


#########################################
#          MAIN                         #
#########################################
reads = pybedtools.BedTool(sys.argv[1])
exons  = pybedtools.BedTool(sys.argv[2])



starts = reads.each(makeStarts).intersect(exons, wo=True).each(convert_to_mRNA)[:]
ends   = reads.each(makeEnds).intersect(exons, wo=True).each(convert_to_mRNA)[:]
while True:
    
    try:
        start = starts.next()
        end   = ends.next()
    except:
        break

    #this still doesn't take into account getting out of sync, so I may loose a lot of reads, I'll take that hit for now, fix later
    if start.name == end.name:
        start.end = end.end
        if start.strand == "-":
            tmp = start.end
            start.end = start.start
            start.start = tmp
        print start,

    
