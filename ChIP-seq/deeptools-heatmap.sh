computeMatrix scale-regions -S XW-GV-H3K4me3.unique_counts_track-1.bw XW-GV-H3K27me3.unique_counts_track-1.bw GV-Shen.unique_counts_track-1.norm.bw -R mm9.TSS.RNA.order.bed -b 10000 -a 10000 -o H2AK119ub.promoter.RNA.Gene.gz --regionBodyLength 10000 --numberOfProcessors 30 --skipZeros --missingDataAsZero


zcat H2AK119ub.promoter.RNA.Gene.gz > H2AK119ub.RNAdecreasing.1


### run in R
x <- read.table("H2AK119ub.RNAdecreasing.1")
x2 <- x[,7:9006]
mean2 <- function(x) {return(mean(x,na.rm=T))}
x3 <- apply(x2,1,mean2)
xout <- x[x3>0.05,]
xout[is.na(xout)] <- "nan"
dim(xout)
write.table(xout,"H2AK119ub.RNAdecreasing.2",col.names=F,row.names=F,sep="\t",quote=F)
###


#@{"upstream":[10000,10000,10000],"downstream":[10000,10000,10000],"body":[10000,10000,10000],"bin size":[10,10,10],"ref point":[null,null,null],"verbose":false,"bin avg type":"mean","missing data as zero":true,"min threshold":null,"max threshold":null,"scale":1,"skip zeros":true,"nan after end":false,"proc number":30,"sort regions":"keep","sort using":"mean","unscaled 5 prime":[0,0,0],"unscaled 3 prime":[0,0,0],"group_labels":["genes"],"group_boundaries":[0,23447],"sample_labels":["XW-GV-H3K4me3.unique_counts_track-1","XW-GV-H3K27me3.unique_counts_track-1","GV-Shen.unique_counts_track-1.norm"],"sample_boundaries":[0,3000,6000,9000]}


plotHeatmap -m H2AK119ub.RNAdecreasing.2.gz -out H2AK119ub.RNAdecreasing.pdf --samplesLabel H3K4me3 H3K27me3 H2AK119ub --regionsLabel "RNAdecreasing" --yMin 0 --colorList "white,blue" --sortRegions no --whatToShow "heatmap and colorbar" --zMax 0.5
