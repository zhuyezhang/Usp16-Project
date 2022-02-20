library(gplots)
library(dplyr)

x <- read.table("Usp16_Promoter.txt")
colnames(x) <- c("chr","start","end","NM","GeneName","strand","H3K4me3","H3K27me3","GV","MII","1C","2C")


table <- x %>% group_by(GeneName) %>% summarise(H3K4me3=mean(H3K4me3),
                                                H3K27me3=mean(H3K27me3),
                                                GV=mean(GV),
                                                MII=mean(MII),
                                                `1C`=mean(`1C`),
                                                `2C`=mean(`2C`))

all <- table
all$MII <- 0.4395839 * all$MII
all$`1C` <- 0.9965277 * all$`1C`
all$`2C` <- 1.0105443 * all$`2C`

library(dplyr)
all <- all %>% mutate(H3K4me3=ifelse(H3K4me3>quantile(H3K4me3,probs = 0.99),quantile(H3K4me3,probs = 0.99) , H3K4me3),
                      H3K27me3=ifelse(H3K27me3>quantile(H3K27me3,probs = 0.99),quantile(H3K27me3,probs = 0.99) , H3K27me3),
                      GV=ifelse(GV>quantile(GV,probs = 0.99),quantile(GV,probs = 0.99) , GV) , 
                      `1C`=ifelse(`1C`>quantile(`1C`,probs = 0.99),quantile(`1C`,probs = 0.99) , `1C`),
                      `2C`=ifelse(`2C`>quantile(`2C`,probs = 0.99),quantile(`2C`,probs = 0.99) , `2C`),
                      `MII`=ifelse(`MII`>quantile(`MII`,probs = 0.99),quantile(`MII`,probs = 0.99) , `MII`))

all <- data.frame(all)
rownames(all) <- all$GeneName
all <- all[,-1]

all <- data.frame(all)

set.seed(2000)
km <- kmeans(all[,1:3],7)
km

library(factoextra)
df <- all[,1:3]
df <- df[sample(1:23284,2000),]

types <- data.frame(Gene=rownames(all),type=km$cluster)
all$Gene <- rownames(all)
allx <- all[,c(7,1,2,3)]
types2 <- merge(allx,types,by="Gene")
library(reshape2)
types3 <- melt(types2,id.vars = c("Gene","type"))
library(ggplot2)
ggplot(types3,aes(x=variable,y=value)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~type) + theme_classic()


expression <- read.table("GSE71434_FPKM_stage.txt",header=T)
expression <- expression[,c(1,4)]
expression <- unique(expression)
expression$X8w_oocyte <- log2(expression$X8w_oocyte + 1)

g1 <- expression[expression$Gene %in% types2$Gene[types2$type==1],]
g2 <- expression[expression$Gene %in% types2$Gene[types2$type==2],]
g3 <- expression[expression$Gene %in% types2$Gene[types2$type==7],]
g4 <- expression[expression$Gene %in% types2$Gene[types2$type==5] | expression$Gene %in% types2$Gene[types2$type==6],]
g5 <- expression[expression$Gene %in% types2$Gene[types2$type==3] | expression$Gene %in% types2$Gene[types2$type==4],]

par(las=2)
boxplot(g1$X8w_oocyte,g2$X8w_oocyte,g3$X8w_oocyte,g4$X8w_oocyte,g5$X8w_oocyte,outline = F,names = c("K27 + ub (n=1,447)","K4 + K27 + ub (n=1,542)","K4 + ub (n=2,507)","K4 only (n=11,110)","None (n=6,678)"))

summary(as.factor(types2$type))


x1 <- types2[types2$Gene %in% types2$Gene[types2$type==1],]
x2 <- types2[types2$Gene %in% types2$Gene[types2$type==2],]
x3 <- types2[types2$Gene %in% types2$Gene[types2$type==7],]
x4 <- types2[types2$Gene %in% types2$Gene[types2$type==5] | types2$Gene %in% types2$Gene[types2$type==6],]
x5 <- types2[types2$Gene %in% types2$Gene[types2$type==3] | types2$Gene %in% types2$Gene[types2$type==4],]

shuffles <- function(x) {
  set.seed(20000)
  x2 <- x[sample(1:dim(x)[1],dim(x)[1]),]
  return(x2)
}

x1 <- shuffles(x1)
x2 <- shuffles(x2)
x3 <- shuffles(x3)
x4 <- shuffles(x4)
x5 <- shuffles(x5)

out <- rbind(x1,x2,x3,x4,x5)
row.names(out) <- out$Gene
out <- out[,c(-1,-5)]
out <- log2(out + 0.1)
out$Gene <- row.names(out)
write.table(out[,c(4,4,1,2,3)],"Cluster2.txt",col.names = T,row.names = F,sep = "\t",quote = F)
