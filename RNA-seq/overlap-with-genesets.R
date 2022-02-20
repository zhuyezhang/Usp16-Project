order1 <- read.table("GV.H3K27me3.H2A.yuan.kmeans.RNA.txt",header = T)
up <- read.table("BW/2cellup.TSS.bed")
down <- read.table("BW/2celldown.TSS.bed")

genes <- c("Zfp352","Zscan4c","Zscan4f","Klf10","Cdk16","Phlda2","Ankrd22","Hspa2")

order1$all <- 1

order1$up <- 0
order1$up[order1$Gene %in% up$V5] <- 1

order1$down <- 0
order1$down[order1$Gene %in% down$V5] <- 1

order1$specific <- 0
order1$specific[order1$Gene %in% genes] <- 1

#write.table(order1[,c(1,2,4,5,6)],"c1.txt",col.names = T,row.names = F,sep="\t",quote = F)

a1 <- order1[1:(4975+6679),]
a2 <- order1[(4975+6679+1):(4975+6679+2959),]
a3 <- order1[(4975+6679+2959+1):(4975+6679+2959+2532),]
a4 <- order1[(4975+6679+2959+2532+1):(4975+6679+2959+2532+6139),]

c1 <- apply(a1[,4:7], 2, sum)
c2 <- apply(a2[,4:7], 2, sum)
c3 <- apply(a3[,4:7], 2, sum)
c4 <- apply(a4[,4:7], 2, sum)

ss <- rbind(c1,c2,c3,c4)

write.table(ss,"c1.txt",row.names = T,col.names = T,sep="\t",quote = F)


ss <- read.table("c1.txt",header = T)
ss <- t(ss)
library(dplyr)
library(reshape2)
library(ggplot2)

ss2 <- melt(ss)
ss3 <- ss2 %>% group_by(Var1) %>% summarise(sums=sum(value))

ss4 <- ss2 %>% left_join(ss3,by="Var1")
ss4$fre <- ss4$value/ss4$sums

ggplot(ss4,aes(x=Var1,y=fre,fill=Var2)) + geom_bar(stat = "identity",width = 0.5) + theme_classic()

#ss4 <- ss4[ss4$Var1 %in% c("all","up","down"),]
#ggplot(ss4,aes(x=Var1,y=fre,fill=Var2)) + geom_bar(stat = "identity",width = 0.5) + theme_classic()
##################################################



order1 <- read.table("Cluster2.txt",header = T)
up <- read.table("BW/2cellup.TSS.bed")
down <- read.table("BW/2celldown.TSS.bed")

genes <- c("Zfp352","Zscan4c","Zscan4f","Klf10","Cdk16","Phlda2","Ankrd22","Hspa2")

order1$all <- 1

order1$up <- 0
order1$up[order1$Gene %in% up$V5] <- 1

order1$down <- 0
order1$down[order1$Gene %in% down$V5] <- 1

order1$specific <- 0
order1$specific[order1$Gene %in% genes] <- 1

#write.table(order1[,c(1,2,4,5,6)],"c2.txt",col.names = T,row.names = F,sep="\t",quote = F)

a1 <- order1[1:(1447),]
a2 <- order1[(1447+1):(1447+1542),]
a3 <- order1[(1447+1542+1):(1447+1542+2507),]
a4 <- order1[(1447+1542+2507+1):(1447+1542+2507+11110),]
a5 <- order1[(1447+1542+2507+11110+1):(1447+1542+2507+11110+6678),]

c1 <- apply(a1[,6:9], 2, sum)
c2 <- apply(a2[,6:9], 2, sum)
c3 <- apply(a3[,6:9], 2, sum)
c4 <- apply(a4[,6:9], 2, sum)
c5 <- apply(a5[,6:9], 2, sum)

ss <- rbind(c1,c2,c3,c4,c5)

write.table(ss,"c2.txt",row.names = T,col.names = T,sep="\t",quote = F)


ss <- read.table("c2.txt",header = T)
ss <- t(ss)
library(dplyr)
library(reshape2)
library(ggplot2)

ss2 <- melt(ss)
ss3 <- ss2 %>% group_by(Var1) %>% summarise(sums=sum(value))

ss4 <- ss2 %>% left_join(ss3,by="Var1")
ss4$fre <- ss4$value/ss4$sums

ggplot(ss4,aes(x=Var1,y=fre,fill=Var2)) + geom_bar(stat = "identity",width = 0.5) + theme_classic()

#ss4 <- ss4[ss4$Var1 %in% c("all","up","down"),]
#ggplot(ss4,aes(x=Var1,y=fre,fill=Var2)) + geom_bar(stat = "identity",width = 0.5) + theme_classic()
