#!/bin/bash

for dir in *.bam
do
dir=${dir/.bam/}
bamToBed -i ${dir}.bam |extend_single ~/tools/ucsc_tools/mm9.chrom.sizes 300|sort -k1,1 -k2,2n > ${dir}.bed
for file in mm9_promoter_2k2k_strand.bed mm9.random.background.bed
do
    file2=${file##*/}
    ./repeats_rpkm.sh ${dir} ${file} ${dir}-${file2/.bed/}.rpkm 
done
done
