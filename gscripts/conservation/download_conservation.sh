#!/bin/bash

#must be manually edited for each species... depending on how many alignments exist.
#only has to be done once though, so that's nice

for PTH in multiz7way phastCons7way phyloP7way
do
   which rsync
   rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/ce10/$PTH/ ./$PTH/;
done