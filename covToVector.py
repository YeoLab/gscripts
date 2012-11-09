#!/usr/bin/python
import sys

vector = {}
prevmRNA = ""
vector = []
for line in open(sys.argv[1]):
    line = line.strip().split()
    curmRNA = line[0]
    if curmRNA != prevmRNA:
        print prevmRNA + "\t" + "\t".join(vector)
        vector = []
        
    vector.append(line[3])
    prevmRNA = curmRNA 

