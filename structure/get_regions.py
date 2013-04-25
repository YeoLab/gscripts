#gets regions from wiggle file
import sys
import bx.bbi.bigwig_file
import pybedtools
from numpy import isnan

def get_genome_sizes(genome_file):
    #input file handle genome_size
    #gets genome sizes indexed as a dict, {chr : int}

    genome = {}

    for line in genome_file:
        line = line.strip().split()
        genome[line[0]] = int(line[1])
    return genome

total_positive_strand = bx.bbi.bigwig_file.BigWigFile(open(sys.argv[1], "rb"))
total_negative_strand = bx.bbi.bigwig_file.BigWigFile(open(sys.argv[2], "rb"))
genome = get_genome_sizes(open(sys.argv[4], 'r'))

tool = pybedtools.BedTool(sys.argv[3])
#count = 0
cur_chrom = None
for line in tool:
#    count += 1
    #assume bed12, will program in standard bed file later...
    #this is for large input sizes (> 1tb)
    #if line.chrom != cur_chrom:
    #    positive_strand = total_positive_strand.get_as_array(line.chrom, 0, genome[line.chrom])
    #    negative_strand = total_negative_strand.get_as_array(line.chrom, 0, genome[line.chrom])
    #    cur_chrom = line.chrom
    
    if len(line.fields) == 12:
        lengths = [ int(x) for x in line.fields[10].split(",") ] 
        start_sites = [ int(x) for x in line.fields[11].split(",") ]
        extra = []

        first_start_site = line.start + start_sites[0]
        first_stop_site = line.start + start_sites[0] + lengths[0]

        second_start_site = line.start + start_sites[1]
        second_stop_site  = line.start + start_sites[1] + lengths[1]
        
        if line.strand == "+":
            print line.chrom, first_start_site, first_stop_site
            first  = total_positive_strand.get_as_array(line.chrom, first_start_site, first_stop_site)
            second = total_positive_strand.get_as_array(line.chrom, second_start_site, second_stop_site)
            #first = positive_strand[first_start_site : first_stop_site ]
            #second = positive_strand[second_start_site : second_stop_site ]
            
        elif line.strand == "-":
            total_negative_strand.get_as_array(line.chrom, first_start_site, first_stop_site)
            total_negative_strand.get_as_array(line.chrom, second_start_site, second_stop_site)
            #first = negative_strand[first_start_site : first_stop_site]    
            #second = negative_strand[second_start_site : second_stop_site]
            
            
        else:
            raise ValueError("strand needs to be either + or - not %s" % (line.strand))

        #first[isnan(first)] = 0
        first = [ str(x) for x in list(first) ]
        #first = [ "0" if x == "0.0" else x for x in first ]

        extra.append(",".join(first))
        
        #second[isnan(second)] = 0
        second = [ str(x) for x in list(second) ]
        #second = [ "0" if x == "0.0" else x for x in second  ]
        extra.append(",".join(second))
        
        #this will get removed once I finish error checking
        line.score = str(int(float(line.score)))
        print str(line).strip() + "\t" + "\t".join(extra)  
    

#        if count > 10000:
#            break
