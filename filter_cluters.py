#!/usr/bin/python
#filters give python clusters to a specific level (bonferoni corrected)

import sys

sig_value = float(sys.argv[2])

#count number of independent tests
count = 0
for line in open(sys.argv[1]):
    count += 1

for line in open(sys.argv[1]):
    line_split = line.split()
    
    
    if not line.startswith("track name") and sig_value / count > float(line_split[4]):
        print line,


