#!/usr/bin/python
'''
Created on Jan 24, 2014

@author: gpratt
'''

import argparse
import os
from string import Template
from gscripts import editing


editing_tempate = os.path.join(editing.editing_dir(), "editing.template")

parser = argparse.ArgumentParser(description="""Takes a bedgraph and an indexed bam file and normalized bedgraph by number of reads in bam file (RPM)
outputs new normalized bedgraph file """)
parser.add_argument("--make_file", help="Make file template to create make file from", 
                    default = editing_tempate, required=False)
parser.add_argument("--bam", help="bam file to create make file from", required=True)
parser.add_argument("--snpEffDb", help="species to create snpEffDb from (hg19 mm9 ect...)", required=True)
parser.add_argument("--snpDb", help="Actual SNP database (from DB SNP)", required=True)
parser.add_argument("--genome", help="fasta genome file to perform snp calling against", required=True)
parser.add_argument("--flipped", help="""if the strands are flipped, not flipped, or not strand specific (specifcy --flipped flipped if it is filipped, both if not strand specifc, and nothing if not flipped)""", default="normal", required=False)
parser.add_argument("--outfile", help="make file used to run rna editing pipeline", required=True)
parser.add_argument("--MinCoverage", help="Minimun Coverage", default="5", required=False)
parser.add_argument("--MinConfidence", help="Minimum Confidence", default="0.995", required=False)
parser.add_argument("--MinEditFrac", help="Mininimum Editing Fraction (G / A+G)", default="0.1", required=False)
parser.add_argument("--pseudoG", help="G Pseudo Counts", required=False, default="5")
parser.add_argument("--pseudoA", help="A Pseudo Counts", required=False, default="5")

args = parser.parse_args()

class MyTemplate(Template):
    delimiter = '&'

with open(args.make_file) as input:
    template = MyTemplate(input.read())

trimmed_name = os.path.splitext(args.bam)[0]
flipped = "S" if args.flipped == "flipped" else "s"

#create simlink to actually run the make file
os.system("ln -s %s %s_stage1.sorted.bam" % (args.bam, trimmed_name))

result = template.substitute(SAMPLE=trimmed_name, 
                             snpEffDB=args.snpEffDb, 
                             SNP_DB=args.snpDb, 
                             S=flipped, 
                             FA=args.genome,
                             MinCoverage=args.MinCoverage,
                             MinConfidence=args.MinConfidence,
                             MinEditFrac=args.MinEditFrac,
                             pseudoG=args.pseudoG,
                             pseudoA=args.pseudoA)

with open(args.outfile, 'w') as out:
    out.write(result)
