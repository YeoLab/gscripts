#!/bin/bash

#"convert to bz2 file, generate a bz2 table, create maf index"
#bzip-table (from seek-bzip2) and maf_build_index.py (from bx-python) must be in PATH

#usage: sh index_conservation.sh chrX.maf.gz

species='ce10'
xx=`basename $1 .gz`
gunzip -c $1 |bzip2 -z > $xx.bz2
bzip-table < $xx.bz2 > $xx.bz2t
maf_build_index.py -s $species $xx.bz2