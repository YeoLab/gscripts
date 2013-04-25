#!/usr/bin/python

#first file is the filtered gene list
#second file is the file to to read and output a filtered list of
import sys

genes = set([])
for gene in open(sys.argv[1]):
    gene = gene.split()[4]
    genes.add(gene)


for line in open(sys.argv[2]):
    gene = line.split()[0]
    if gene in genes:
        print line,
