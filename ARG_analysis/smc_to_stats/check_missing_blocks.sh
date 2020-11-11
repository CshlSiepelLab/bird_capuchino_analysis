#!/bin/bash

# simple bash script for checking for missing ARG blocks after window stats have been computed

window_info_file=./infoTables/scaffold-window-coverage.txt
scaffold_length_file=  # <-- SET FILE WITH LENGTH FOR EACH SCAFFOLD / CHROMOSOME
# effective min scaffold length is 110000, because we need at least one 20 kb window after trimming
min_scaffold_len=110000
# left trim should be accurate and right trim should be <= trim + 1/2 window length
trim_left=50000
max_trim_right=60000

# filter scaffolds based on length
cat $scaffold_length_file | awk -v min_scaffold_len=$min_scaffold_len '{if($2>=min_scaffold_len){print;}}' | sort -k1,1 > temp_file.txt

sort -k1,1 $window_info_file | join - temp_file.txt -a 1 -a 2 > temp_file2.txt

mv temp_file2.txt temp_file.txt

cat temp_file.txt | awk -v trim_left=$trim_left -v max_trim_right=$max_trim_right '{if(NF!=4 || $2 != trim_left || $4-$3 > max_trim_right) {print;}}'
rm temp_file.txt
