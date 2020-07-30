#!/bin/bash

# bash script to compute scaffolds that do not contain F_ST peaks
#  and split them into blocks of given size that act as a control
#  set for all analyses of Fst peaks. These should be long enough
#  s.t. the great majority of peaks are shorter.

FSTfile=infoTables/FST-table-base.txt
coverageFile=infoTables/scaffold-window-coverage.txt
outFile=infoTables/control-scaffolds.txt
outFile_blocks=infoTables/control-blocks.txt
blockLen=500000      #  length of blocks you want to act as controls for entire Fst peaks (our set of 1,376 blocks)


cut -f2 $FSTfile | tail -n+2 | sort | uniq | sed -r 's|^.*$|Contig& FST|' > $outFile.tmp
# echo $fst_scaffolds
 
cat $coverageFile | grep -v "^#" | sed -r 's|^Contig([0-9]+)[[:space:]]+|Contig\1 \1 |' | sort -k1,1 | join - $outFile.tmp -a1 | grep -v "FST$" | sort -k2,2n | cut -d" " -f1,3- > $outFile

num_fst_scaffolds=`grep -c . $outFile.tmp`
num_control_scaffolds=`grep -c . $outFile`
num_scaffolds=`grep -vc "^#" $coverageFile`

echo "Found $num_fst_scaffolds F_ST scaffolds and $num_control_scaffolds control scaffolds ($num_scaffolds total)"
rm $outFile.tmp

# find long scaffolds and split into blocks

cat $outFile | awk -v blockLen=$blockLen \
  '{ \
     blockStart = $2; \
     blockEnd   = blockStart + blockLen; \
     while(blockEnd <= $3) { \
       print $1,blockStart,blockEnd; \
       blockStart = blockEnd; \
       blockEnd   = blockStart + blockLen; \
     } \
  }' > $outFile_blocks

num_control_blocks=`grep -c . $outFile_blocks`

echo "Found $num_control_blocks control blocks of size $blockLen in these scaffolds"
