############################

library(ggplot2)
require(gtable)
library(grid)
library(ggplot2)
library(dplyr)
library(reshape2)


expression_table <- read.table("genes.fpkm_table",header = T)

expression_table <- expression_table[-(grep("^Mir",expression_table[,"tracking_id"])),]
expression_table <- expression_table[-(grep("^Sno",expression_table[,"tracking_id"])),]
#expression_table <- expression_table[-(grep("^ERCC",expression_table[,"tracking_id"])),]

fpkm_1 <-c()
fpkm_5 <-c()
fpkm_10 <-c()
for (i in 2:dim(expression_table)[2]){
  fpkm_1 <- cbind(fpkm_1,sum(expression_table[,i]>1))
  fpkm_5 <- cbind(fpkm_5,sum(expression_table[,i]>5))
  fpkm_10 <- cbind(fpkm_10,sum(expression_table[,i]>10))
}

fpkm_genes <- rbind(fpkm_1,fpkm_5,fpkm_10)
colnames(fpkm_genes) <- colnames(expression_table)[2:dim(expression_table)[2]]
rownames(fpkm_genes) <- c("FPKM>1","FPKM>5","FPKM>10")
fpkm_genes <- data.frame(fpkm_genes)

outputfile_genenumber <- file.path("Expression_genenumber_byFPKM.csv")
#write.csv(fpkm_genes,outputfile_genenumber,quote = T,row.names = T)

rownames(expression_table) <- expression_table$tracking_id
expression_table <- expression_table[,c(1,8,10,2,3,12,13,5,6)]

expression_table <- expression_table[,-1]

number_greater_than1 <- function(vector_a) {
  return(sum(vector_a>1))
}

expression_table$number_greater_than1 <- apply(expression_table,1,number_greater_than1)
expression_table <- expression_table[expression_table[,dim(expression_table)[2]]>=1,]
expression_table <- expression_table[,-(dim(expression_table)[2])]
expression_table[expression_table<1] = 1

write.table(expression_table,"Usp16_genes_FPKM.txt",col.names = T,row.names = T,sep = "\t",quote = F)

library(corrplot)
col1 <- colorRampPalette(c("green","yellow","blue","white","red"))

corrplot(cor(expression_table,method = "spearman"),col=col1(51),method = "color",cl.lim = c(0,1),is.corr=FALSE,tl.col = "black",addCoef.col = "black")

expression_table$WT1C <- 0.5 * (expression_table$WT.1cell.r1_0 + expression_table$WT.1cell.r3_0)
expression_table$Usp16KO1C <- 0.5 * (expression_table$Usp16.cKO.1cell.r1_0 + expression_table$Usp16.cKO.1cell.r2_0)
expression_table$WT2C <- 0.5 * (expression_table$WT.2cell.r2_0 + expression_table$WT.2cell.r3_0)
expression_table$Usp16KO2C <- 0.5 * (expression_table$Usp16.cKO.2cell.r1_0 + expression_table$Usp16.cKO.2cell.r2_0)

expression_table <- expression_table[,c(9,10,11,12)]

oup <- expression_table[expression_table$Usp16KO1C > 3 * expression_table$WT1C,]
odown <- expression_table[expression_table$WT1C > 3 * expression_table$Usp16KO1C,]

onochanged <-expression_table[!rownames(expression_table) %in% rownames(oup) & !rownames(expression_table) %in% rownames(odown),]
oup$type <- "Up"
odown$type <- "Down"
onochanged$type <- "Nochange"

plot_genes <- rbind(oup,odown,onochanged)
library(ggplot2)
plot_genes$type <- factor(plot_genes$type,levels=c("Down","Nochange","Up"))
ggplot(plot_genes,aes(x=log2(WT1C),y=log2(Usp16KO1C),col=type)) + geom_point() + scale_color_manual(values = c("blue","grey","red")) +
  theme_classic() + theme(legend.position = 'none') + scale_y_continuous(limits = c(0,16)) + scale_x_continuous(limits = c(0,16))

tup <- expression_table[expression_table$Usp16KO2C > 3 * expression_table$WT2C,]
tdown <- expression_table[expression_table$WT2C > 3 * expression_table$Usp16KO2C,]

tnochanged <-expression_table[!rownames(expression_table) %in% rownames(tup) & !rownames(expression_table) %in% rownames(tdown),]
tup$type <- "Up"
tdown$type <- "Down"
tnochanged$type <- "Nochange"

plot_genes <- rbind(tup,tdown,tnochanged)
library(ggplot2)
plot_genes$type <- factor(plot_genes$type,levels=c("Down","Nochange","Up"))
ggplot(plot_genes,aes(x=log2(WT2C),y=log2(Usp16KO2C),col=type)) + geom_point() + scale_color_manual(values = c("blue","grey","red")) +
  theme_classic() + theme(legend.position = 'none') + scale_y_continuous(limits = c(0,16)) + scale_x_continuous(limits = c(0,16))


zga2 <- expression_table[expression_table$WT2C > 3 * expression_table$WT1C & expression_table$WT2C > 1,]
library(VennDiagram)
T1<-venn.diagram(x =list(ZGA = rownames(zga2),`2CellDown` = rownames(tdown)), filename = NULL, height = 600, width= 600, resolution =150, imagetype="png", fill =c("green","red"),alpha=0.5,cex=2.5,cat.cex=2.5,cat.pos=c(0, 0))
grid.draw(T1)

write.table(oup,"1Cell-up.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(tup,"2Cell-up.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(odown,"1Cell-down.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(tdown,"2Cell-down.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(zga2,"ZGA.defined.txt",col.names = T,row.names = T,sep = "\t",quote = F)

phyper(603,603+1337,19000,603+3136,log.p = T,lower.tail = F)/log(10)
