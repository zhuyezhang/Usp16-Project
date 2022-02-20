library(chromstaR)

files <- list.files('./',"Shen.unique.bam$")

mm9_chrominfo <- read.table("mm9.chrom.sizes",header = F)
colnames(mm9_chrominfo) <- c("chromosome","length")

for (file in files){
file2 <- gsub(".bam","",file)
binned.data <- binReads(file,assembly = mm9_chrominfo,binsizes = 1000,stepsizes = 200)
model <- callPeaksUnivariate(binned.data = binned.data,verbosity = 0)
exportCounts(model,filename = file2)
exportPeaks(model,filename = file2)
}

