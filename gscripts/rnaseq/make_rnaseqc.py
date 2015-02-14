#!/usr/bin/python
#makes and outputs to stdout a formatted file to run rnaseq-qc on
import sys
import os

for sample in sys.argv[1:]:
    print "\t".join(["Sample ID", "Bam File",    "Notes"])
    print "\t".join([os.path.basename(sample).split(".")[0], os.path.abspath(sample), "foo"])
