__author__ = 'lovci'


import argparse

parser = argparse.ArgumentParser(description='Build an index suitable for gsnap, /'
                                             'http://research-pub.gene.com/gmap/src/README')

parser.add_argument('--index_name',
                   help='a name for the index')

parser.add_argument('--files', action='append',
                   help='fasta files to be put in the index')

opts = parser.parse_args()
import subprocess

exit(subprocess.call(['gmap_build', '-d', opts.index_name, '-k', '15', " ".join(opts.files)]))