GENOME=/nethome/Shared/annotation/iGenome/Mus_musculus/UCSC/mm9/Sequence/Bowtie2Index/genome

# obtain first 50bp of reads
for file in `ls *R1*.fastq.gz`
do
{
trim-gz-c $file ${file%%_*}_1.fastq.gz 0 50
}&
done
wait


# mapping
for file in `ls *_1.fastq.gz`
do
bowtie2 -p 12 -x $GENOME -U $file -S ${file/_1.fastq.gz/}.sam 2> ${file/_1.fastq.gz/}.align.log
grep -v "XS:i:" ${file/_1.fastq.gz/}.sam > ${file/_1.fastq.gz/}.unique.sam
rm ${file/_1.fastq.gz/}.sam
done

# peak calling
for file in `ls *.sam`
do
file2=${file#*/}
file3=${file2/.sam/}
macs2 callpeak -f SAM -q 0.05 --nomodel --nolambda --broad --extsize 300 -B --SPMR -g mm --outdir $file3 -n $file3 -t $file
done

# peak distribution
for file in *bed 
do
ceasBW -w ${file/.bed/.bed} -b $file -g ~/tools/ucsc_tools/mm9.refGene -l ~/tools/ucsc_tools/mm9.chrom.sizes
done