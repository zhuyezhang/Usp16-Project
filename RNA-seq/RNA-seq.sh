#!/bin/bash

GTF=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf
GENOME_FA=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome.fa
GENOME=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

for file in `ls *R1*.fastq`
do 
{
trim-c $file ${file%%_*}_1.fastq 0 50
} &
done
wait

for file in `ls *_1.fastq`
do
outdir=${file%_*}_tophat
tophat -o ${outdir} -p 20 -G ${GTF} ${GENOME} ${file} #--no-novel-juncs 
done

for filedir in `ls -F | grep 'tophat/$'`
do
{
bash ~/tools/mycode/filter_unique_reads_from_tophat.sh ${filedir}accepted_hits.bam ${filedir}accepted_hits_unique
}&
done
wait


for filedir in `ls -F | grep 'tophat/$'`
do
outdir=${filedir%%/}_cuffquant
cuffquant -o ${outdir} -p 16 -u $GTF ${filedir}accepted_hits_unique.bam
done

files2=""
label=""
for file in `ls -F | grep '_cuffquant/$'`
do
label=${label}${file%%_*}","
files2=${files2}${file}"abundances.cxb "
done

label=${label%,*}

cuffnorm -p 20 --labels ${label} --library-norm-method=classic-fpkm -o cuffnorm_table $GTF ${files2}

mkdir tophat cuffquant
mv *_tophat/ tophat
mv *_cuffquant/ cuffquant
