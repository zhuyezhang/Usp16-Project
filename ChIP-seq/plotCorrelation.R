table <- read.table("PearsonCorr_results.tab")

table <- table[c(5,6,7,8,1,4,2,3),c(5,6,7,8,1,4,2,3)]

library(corrplot)

s <- c("GV-1","GV-2","MII-1","MII-2","1C-1","1C-2","2C-1","2C-2")
colnames(table) <- s
rownames(table) <- s

col1 <- colorRampPalette(c("green","yellow","blue","white","red"))

corrplot(as.matrix(table),col=col1(51),method = "color",order = "original",cl.lim = c(0,1),is.corr=FALSE,tl.col = "black",addCoef.col = "black")




table <- read.table("PromoterPearsonCorr_results.tab")

table <- table[c(2,6,7,8,1,5,3,4),c(2,6,7,8,1,5,3,4)]

library(corrplot)

s <- c("GV-1","GV-2","MII-1","MII-2","1C-1","1C-2","2C-1","2C-2")
colnames(table) <- s
rownames(table) <- s

col1 <- colorRampPalette(c("green","yellow","blue","white","red"))

corrplot(as.matrix(table),col=col1(51),method = "color",order = "original",cl.lim = c(0,1),is.corr=FALSE,tl.col = "black",addCoef.col = "black")



table <- read.table("TSSPearsonCorr_results.tab")

table <- table[c(1,8,2,3,4,7,5,6),c(1,8,2,3,4,7,5,6)]

library(corrplot)

s <- c("GV-1","GV-2","MII-1","MII-2","1C-1","1C-2","2C-1","2C-2")
colnames(table) <- s
rownames(table) <- s

col1 <- colorRampPalette(c("green","yellow","blue","white","red"))

corrplot(as.matrix(table),col=col1(51),method = "color",order = "original",cl.lim = c(0,1),is.corr=FALSE,tl.col = "black",addCoef.col = "black")
