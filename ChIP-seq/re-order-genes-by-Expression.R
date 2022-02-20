genes <- read.table("GSE71434_FPKM_stage.txt",header=T)
genes <- genes[order(genes$X8w_oocyte,decreasing = T),]

tss <- read.table("mm9.TSS.bed")

genes <- data.frame(GeneName=unique(genes$Gene))
colnames(tss) <- c("chr","start","end","NM","GeneName","strand")

library(dplyr)
genes2 <- genes %>% left_join(tss,by="GeneName")
genes2 <- na.omit(genes2)

genes2 <- genes2[,c(2,3,4,5,1,6)]
write.table(genes2,"mm9.TSS.RNA.order.bed",sep="\t",col.names = F,row.names = F,quote = F)

