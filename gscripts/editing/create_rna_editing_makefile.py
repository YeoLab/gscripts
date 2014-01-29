#!/usr/bin/python
'''
Created on Jan 24, 2014

@author: gpratt
'''

import argparse
import os
from string import Template

parser = argparse.ArgumentParser(description="""Takes a bedgraph and an indexed bam file and normalized bedgraph by number of reads in bam file (RPM)
outputs new normalized bedgraph file """)
parser.add_argument("--make_file", help="Make file template to create make file from", 
                    default=os.path.join(os.path.dirname(__file__), "template.mk"), required=False)
parser.add_argument("--bam", help="bam file to create make file from", required=True)
parser.add_argument("--snpEffDb", help="species to create snpEffDb from (hg19 mm9 ect...)", required=True)
parser.add_argument("--snpDb", help="Actual SNP database (from DB SNP)", required=True)
parser.add_argument("--genome", help="fasta genome file to perform snp calling against", required=True)
parser.add_argument("--flipped", help="""if the strands are flipped, not flipped, or not strand specific (specifcy --flipped flipped if it is filipped, 
                                         both if not strand specifc, and nothing if not flipped)""", default="normal", required=False)
parser.add_argument("--outfile", help="make file used to run rna editing pipeline", required=True)

args = parser.parse_args()

class MyTemplate(Template):
    delimiter = '&'

with open(args.make_file) as input:
    template = MyTemplate(input.read())

flipped = "S" if args.flipped == "flipped" else "s"

result = template.substitute(sample=os.path.splitext(args.bam)[0], snpEffDb=args.snpEffDb, snpDb=args.snpDb, flipped=flipped, FA=args.genome, )

with open(args.outfile, 'w') as out:
    out.write(result)
