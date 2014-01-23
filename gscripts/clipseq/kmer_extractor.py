#!/usr/bin/python
    
import os
import cPickle as pickle
import argparse


def motif_to_str(motif):
    return "\t".join(map(str, [motif.freq1, motif.freq2, motif.delta]))

parser = argparse.ArgumentParser(description="""Take kmers from clip analysis pickle file and outputs them to tsv files
one for each region and kmer size""")
parser.add_argument("--input", help="Input pickle file from clip_analysis", required=True)
parser.add_argument("--output", help="output file to print results to", required=True)
args = parser.parse_args()

results = pickle.load(open(args.input))
for region_name, region in results['kmer_results'].items():
    for kmer_name, kmer in region.items():
        out_name, ext = os.path.splitext(args.output)
        with open("%s_%s_%s%s" % (out_name, region_name, str(kmer_name), ext), 'w') as outfile:
            for motif, motif_object in kmer.items():
                outfile.write(motif + "\t" + motif_to_str(motif_object) + "\n")
