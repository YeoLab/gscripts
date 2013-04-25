#!/usr/bin/python

#takes in a vector, outputs any location where the peaks are 25 fold higher than the median of the vector
#Eventually want to modify this to accept a) bed formatted files b) gene body files that let me limit searches
#arg1 = vector
#arg2 = gene body information
import sys

#get gene bodies 
geneBodies = {}
for line in open(sys.argv[2]):
    line = line.strip().split()
    geneBodies[line[0]] = map(int, line[1:])
    
for line in open(sys.argv[1]):
    line = line.strip().split()
    name = line[0]
    if name in geneBodies:
        vector = map(int, line[1:])
        geneBody = vector[geneBodies[name][0]:geneBodies[name][1]]
        #geneBody.sort()
        median = geneBody[len(geneBody) /2]
        if median == 0: 
            median = 1
        for location, value in enumerate(vector):
            
            if (value / median) >= 25:
                print name + "\t" + str(location) + "\t" + str(location + 1)
                #print geneBody[-1] / median         
            
