#R 3.5

library(DESeq2)

promoter <- read.table("PROMOTER.txt")

promoter$id <- paste(promoter$V1,promoter$V2,promoter$V3,promoter$V4,promoter$V5,promoter$V6,sep=":")
promoter2 <- promoter[,c(5,10,14,18,22,26,30,34,38)]

factors2 <- 1/ c(1,0.4395839,0.9965277,1.0105443,1,0.4395839,0.9965277,1.0105443)

colnames(promoter2) <- c("Gene","GV1","MII1","X1C1","X2C1","GV2","MII2","X1C2","X2C2")


library(dplyr)
promoter2 <- promoter2 %>% group_by(Gene) %>% summarise(GV1=round(mean(GV1)),MII1=round(mean(MII1)),X1C1=round(mean(X1C1)),X2C1=round(mean(X2C1)),
                                                        GV2=round(mean(GV2)),MII2=round(mean(MII2)),X1C2=round(mean(X1C2)),X2C2=round(mean(X2C2)))


promoter2 <- data.frame(promoter2)

rownames(promoter2) <- promoter2$Gene

promoter2 <- promoter2[,-1]

data2 <- rbind(promoter2)

coldata <- data.frame(X=c("GV1","MII1","X1C1","X2C1","GV2","MII2","X1C2","X2C2"),condition=c("GV","MII","X1C","X2C","GV","MII","X1C","X2C"),type=c("GV","MII","X1C","X2C","GV","MII","X1C","X2C"))
rownames(coldata) <- coldata$X
coldata <- coldata[,-1]

coldata$condition <- relevel(coldata$condition, ref = "GV")


dds <- DESeqDataSetFromMatrix(countData = data2,
                              colData = coldata,
                              design = ~ condition)

#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

dds2 <- estimateSizeFactors(dds)
#dds2$sizeFactor
normCounts = as.data.frame(counts(dds2, normalized=TRUE))


formersf <- c(15406326,6088183,13834018,18580919,41091613,1867108,24494997,24356979) / 10000000


sizeFactors(dds2)  <-  formersf * factors2

normCounts = as.data.frame(counts(dds2, normalized=TRUE))
normCounts$Gene <- rownames(normCounts)

normCounts <- normCounts[,c(9,1,2,3,4,5,6,7,8)]
write.table(normCounts,"normCounts.txt",col.names = T,row.names = F,sep = "\t",quote=F)


dds_2_lrt <- DESeq(dds2)
resultsNames(dds_2_lrt)


res_2_lrt <- results(dds_2_lrt,name = "condition_X2C_vs_GV")

plotMA(res_2_lrt, ylim=c(-2,2))


res_2_lrt <- data.frame(res_2_lrt)
res_2_lrt$Gene <- rownames(res_2_lrt)

res_2_lrt <- merge(normCounts,res_2_lrt,by="Gene")

res_2_lrt$GV <- 0.5 * (res_2_lrt$GV1 + res_2_lrt$GV2)
res_2_lrt$X2C <- 0.5 * (res_2_lrt$X2C1 + res_2_lrt$X2C2)

down <- data.frame(res_2_lrt[res_2_lrt$log2FoldChange < -0.5 & res_2_lrt$pvalue < 0.1 & res_2_lrt$GV > 60, ])
up <- data.frame(res_2_lrt[res_2_lrt$log2FoldChange > 0.5 & res_2_lrt$pvalue < 0.1 & res_2_lrt$X2C > 60, ])

write.table(down,"2C_vs_GV_down.txt",col.names = T,row.names = F,sep = "\t",quote=F)
write.table(up,"2C_vs_GV_up.txt",col.names = T,row.names = F,sep = "\t",quote=F)





res_2_lrt <- results(dds_2_lrt,name = "condition_MII_vs_GV")

plotMA(res_2_lrt, ylim=c(-2,2))

res_2_lrt <- data.frame(res_2_lrt)
res_2_lrt$Gene <- rownames(res_2_lrt)

res_2_lrt <- merge(normCounts,res_2_lrt,by="Gene")

res_2_lrt$GV <- 0.5 * (res_2_lrt$GV1 + res_2_lrt$GV2)
res_2_lrt$MII <- 0.5 * (res_2_lrt$MII1 + res_2_lrt$MII2)

down <- data.frame(res_2_lrt[res_2_lrt$log2FoldChange < -0.5 & res_2_lrt$pvalue < 0.1 & res_2_lrt$GV > 60, ])
up <- data.frame(res_2_lrt[res_2_lrt$log2FoldChange > 0.5 & res_2_lrt$pvalue < 0.1 & res_2_lrt$MII > 60, ])

write.table(down,"MII_vs_GV_down.txt",col.names = T,row.names = F,sep = "\t",quote=F)
write.table(up,"MII_vs_GV_up.txt",col.names = T,row.names = F,sep = "\t",quote=F)



###################################################################


genes <- read.table("GSE71434_FPKM_stage.txt",header=T)
genes <- unique(genes)
genes <- genes[,c("Gene","X8w_oocyte","MII_oocyte","X2cell_early","X2cell_late")]
tcdown <- read.table("2C_vs_GV_down.txt",header = T)
genes_tcdown <- genes[genes$Gene %in% tcdown$Gene,]

miidown <- read.table("MII_vs_GV_down.txt",header = T)
genes_miidown <- genes[genes$Gene %in% miidown$Gene,]

#other <- genes[!genes$Gene %in% tcdown$Gene & !genes$Gene %in% miidown$Gene,]
set.seed(20000)
randomother <- genes[sample(1:24028,2500),]

genes_tcdown$type <- "2Cdown"
genes_miidown$type <- "MIIdown"
randomother$type <- "Other"

library(reshape2)
genes_komiiup2 <- rbind(genes_tcdown,genes_miidown,randomother)
gene2 <- melt(genes_komiiup2)
library(ggplot2)
ggplot(gene2,aes(x=type,y=log2(value+1),fill=variable)) + geom_boxplot(outlier.shape=NA) + coord_cartesian(ylim=c(0,8.5))  + xlab("") + ylab("log2(FPKM+1)") + theme_classic()

t.test(log2(genes_tcdown$X8w_oocyte+1),log2(genes_tcdown$X2cell_late+1),paired = T)
t.test(log2(genes_miidown$X8w_oocyte+1),log2(genes_miidown$X2cell_late+1),paired = T)
t.test(log2(randomother$X8w_oocyte+1),log2(randomother$X2cell_late+1),paired = T)



summary(genes_tcdown$X8w_oocyte)

getmean <- function(x){
  g1 <- x[x$X8w_oocyte<=0.07176,]
  g2 <- x[x$X8w_oocyte>0.07176 & x$X8w_oocyte<=0.67230,]
  g3 <- x[x$X8w_oocyte>0.67230,]
  r1 <- g1[sample(1:dim(g1)[1],250),]
  r2 <- g2[sample(1:dim(g2)[1],125),]
  r3 <- g3[sample(1:dim(g3)[1],125),]
  r <- rbind(r1,r2,r3)
  return(r)
}

num <- 0
for (x in c(1:10000)) {
  x2 <- getmean(genes)
  av2 <- mean(log2(x2$X2cell_late+1))
  if (av2 > 1.662) {
    num <- num + 1
  }
}

num
x2$type <- "Random"

library(reshape2)
genes_komiiup2 <- rbind(genes_tcdown,x2)
gene2 <- melt(genes_komiiup2)
library(ggplot2)
ggplot(gene2,aes(x=type,y=log2(value+1),fill=variable)) + geom_boxplot(outlier.shape=NA) + coord_cartesian(ylim=c(0,6.5))  + xlab("") + ylab("log2(FPKM+1)") + theme_classic()


###################################################################
x1 <- read.table("2C_vs_GV_down.txt",header = T)

x2 <- read.table("MII_vs_GV_down.txt",header = T)

library(VennDiagram)
T1<-venn.diagram(x =list(`2Celldown` = x1$Gene,`MIIdown` = x2$Gene), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("green","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0, 0))
grid.draw(T1)

