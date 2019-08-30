#!/bin/bash

# Split gen files into individual SNPs
# SYNTAX -l # lines, -d numeric suffixes, -a # number of digits in the suffix 
for file in *gen; do split -l 1 -d -a 4  ${file} ./split/split_${file}; done

# Replace space for new line in each file (transpose)
# tr replace all space ' ' for new line '\n' 
for file in ./split/*; do tr ' ' '\n' < ${file} > ${file}_trans; done

# Drop every third digit using awk
for file in ./split/*trans; do awk 'NR % 3 != 0' ${file} > ${file}2; done

# Merge all the SNP files together
cp split_final_full_snps_1.gen0000_trans2 final_snps.txt
for file in ./split/*trans2; do echo "${file}"; paste final_snps.txt ${file} >> final_snps.txt; done
