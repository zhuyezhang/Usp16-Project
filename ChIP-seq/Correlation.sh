multiBamSummary bins --bamfiles *.bam -out bins.results.npz --extendReads 300 --ignoreDuplicates -p 30
plotCorrelation -in bins.results.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o heatmap_PearsonCorr_results.png --outFileCorMatrix PearsonCorr_results.tab --skipZeros --removeOutliers 


multiBamSummary BED-file --BED mm9_promoter_2k2k_strand.bed --bamfiles *.bam -out Promoter.results.npz --extendReads 300 --ignoreDuplicates -p 30
plotCorrelation -in Promoter.results.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o Promoterheatmap_PearsonCorr_results.png --outFileCorMatrix PromoterPearsonCorr_results.tab --skipZeros --removeOutliers 


multiBamSummary BED-file --BED mm9.TSS.bed --bamfiles *.bam -out TSS.results.npz --extendReads 300 --ignoreDuplicates -p 30
plotCorrelation -in TSS.results.npz --corMethod pearson --skipZeros --plotTitle "Pearson Correlation of Read Counts" --whatToPlot heatmap --colorMap RdYlBu --plotNumbers -o TSSheatmap_PearsonCorr_results.png --outFileCorMatrix TSSPearsonCorr_results.tab --skipZeros --removeOutliers 


