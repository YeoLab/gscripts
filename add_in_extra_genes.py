#!/usr/bin/python

import sys

#argv[1] bed file to generate coverage from

#argv[2] genome file

total_genes = set([])
genes_in_set = set([])
for line in open(sys.argv[1]):
    print line,
    line = line.split("\t")
    genes_in_set.add(line[0])
    
for line in open(sys.argv[2]):
    total_genes.add(line.split("\t")[0])

leftover_genes = total_genes - genes_in_set

for leftover_gene in leftover_genes:
    print "%s\t%s\t%s" % (leftover_gene,"1","2") 

                    
