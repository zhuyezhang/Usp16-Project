if [ $# -lt 3 ];then
echo "Need 2 parameters! <IP_sample> <peaklist> <outfile>"
exit
fi

IP_sample=$1
outfile=$3
IP_bam=$(echo "${IP_sample}")
IP_mappable_read_count=$(samtools view -F 0x0004 ${IP_bam}.bam | wc -l)


peakfile=$2

cut -f 1,2,3 ${peakfile}|coverageBed -b ${IP_sample}.bed -a - | awk -v OFS='\t' '{print $1,$2,$3,$4}' > $outfile
