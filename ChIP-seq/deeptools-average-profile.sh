computeMatrix reference-point -S GV-Shen.unique_counts_track-1.norm.bw MII-Shen.unique_counts_track-1.norm.bw 1C-Shen.unique_counts_track-1.norm.bw 2C-Shen.unique_counts_track-1.norm.bw -R 2C_vs_GV_down.bed -b 10000 -a 10000 -o H2AK119ub.promoter.2Cdown.Gene.gz --numberOfProcessors 30 --referencePoint center

plotProfile -m H2AK119ub.promoter.2Cdown.Gene.gz -out H2AK119ub.promoter.2Cdown.Gene.pdf --samplesLabel GV MII 1Cell 2Cell --regionsLabel "2CellDown (n=329)" --yMin 0 --plotHeight 5 --plotWidth 5 --perGroup



computeMatrix reference-point -S GV-Shen.unique_counts_track-1.norm.bw MII-Shen.unique_counts_track-1.norm.bw 1C-Shen.unique_counts_track-1.norm.bw 2C-Shen.unique_counts_track-1.norm.bw -R MII_vs_GV_down.bed -b 10000 -a 10000 -o H2AK119ub.promoter.MIIdown.Gene.gz --numberOfProcessors 30 --referencePoint center

plotProfile -m H2AK119ub.promoter.MIIdown.Gene.gz -out H2AK119ub.promoter.MIIdown.Gene.pdf --samplesLabel GV MII 1Cell 2Cell --regionsLabel "MIIDown (n=2288)" --yMin 0 --plotHeight 5 --plotWidth 5 --perGroup




computeMatrix scale-regions -S GV-Shen.unique_counts_track-1.norm.bw 2C-Shen.unique_counts_track-1.norm.bw -R ZGA.defined.TSS.bed Random.defined.bed -b 10000 -a 10000 -o H2AK119ub.promoter.ZGAdefined.Gene.gz --regionBodyLength 10000 --numberOfProcessors 30

plotProfile -m H2AK119ub.promoter.ZGAdefined.Gene.gz -out H2AK119ub.promoter.ZGAdefined.Gene.pdf --samplesLabel GV 2Cell --regionsLabel "ZGA defined" "Random" --yMin 0 --plotHeight 5 --plotWidth 5 --perGroup



computeMatrix scale-regions -S GV-Shen.unique_counts_track-1.norm.bw 2C-Shen.unique_counts_track-1.norm.bw -R Maternal.bed Minor.bed Major.bed -b 10000 -a 10000 -o H2AK119ub.promoter.DBTME.Gene.gz --regionBodyLength 10000 --numberOfProcessors 30

plotProfile -m H2AK119ub.promoter.DBTME.Gene.gz -out H2AK119ub.promoter.DBTME.Gene.pdf --samplesLabel GV 2Cell --regionsLabel "Maternal" "Minor" "Major" --yMin 0 --plotHeight 5 --plotWidth 5 --perGroup
