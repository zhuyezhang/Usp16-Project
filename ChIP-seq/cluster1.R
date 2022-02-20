# R 4.0.5

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

library(factoextra)
df <- all[,1:2]
df <- df[sample(1:23284,2000),]

set.seed(1000)
km <- kmeans(all[,1:2],5)
km
types <- data.frame(Gene=rownames(all),type=km$cluster)
all$Gene <- rownames(all)
allx <- all[,c(7,1,2)]
types2 <- merge(allx,types,by="Gene")
library(reshape2)
types3 <- melt(types2,id.vars = c("Gene","type"))
library(ggplot2)
#ggplot(types3,aes(x=variable,y=value)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~type)

types2 <- merge(all,types,by="Gene")
types2$type[types2$type==1] <- "None"
types2$type[types2$type==2] <- "Bivalent"
types2$type[types2$type==3] <- "H3K4me3-low"
types2$type[types2$type==4] <- "H3K4me3-high"
types2$type[types2$type==5] <- "H3K27me3"

c1 <- types2[types2$type=="H3K4me3-high",]
c2 <- types2[types2$type=="H3K4me3-low",]
c3 <- types2[types2$type=="Bivalent",]
c4 <- types2[types2$type=="H3K27me3",]
c5 <- types2[types2$type=="None",]

out <- rbind(c1,c2,c3,c4,c5)

out2 <- out[,c(1,2,3,8)]
types3 <- melt(out2,id.vars = c("Gene","type"))
library(ggplot2)
ggplot(types3,aes(x=variable,y=value)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~type) + xlab("") + theme_classic()


all2 <- out[,c(1,1,2,3,4,5,6,7)]

write.table(all2,"GV.H3K27me3.H2A.yuan.kmeans.txt",quote=F,sep="\t",col.names = T,row.names = F)


out2 <- out[,c(1,2,3,4,5,6,7,8)]
types3 <- melt(out2,id.vars = c("Gene","type"))
types3$type[types3$type=="H3K4me3-high"] <- "H3K4me3"
types3$type[types3$type=="H3K4me3-low"] <- "H3K4me3"
library(ggplot2)
types3$type <- factor(types3$type,levels=c("H3K4me3","H3K27me3","Bivalent","None"))
ggplot(types3,aes(x=variable,y=value)) + geom_boxplot(outlier.shape = NA) + facet_wrap(~type) + theme_classic() + xlab("")



##############################################################################################

expression <- read.table("GSE71434_FPKM_stage.txt",header=T)
expression <- expression[,c(1,4)]
expression <- unique(expression)
expression$X8w_oocyte <- log2(expression$X8w_oocyte + 0.1)

out2 <- all2 %>% left_join(expression,by="Gene")
write.table(out2[,c(1,2,9)],"GV.H3K27me3.H2A.yuan.kmeans.RNA.txt",quote=F,sep="\t",col.names = T,row.names = F)

