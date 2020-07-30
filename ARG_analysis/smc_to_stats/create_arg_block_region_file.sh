#!/bin/bash
# script reads smc.gz files and extracts for each one: the ARG file name, scaffold name, start coordinate, end coordinate
# prints table to output file 
smcDir=  # <-- SET SOURCE DIR CONTAINING SMC FILES
outFile=./infoTables/ARGblock-coordinates.txt

rm -f $outFile
files=`ls $smcDir/*.smc.gz`

for f in $files; do
  line=$f"\t"`zcat $f | head -n2 | tail -n1 | tr -s [[:space:]] "\t" | cut -f2,3,4`
  echo -e $line >> $outFile
done
