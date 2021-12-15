setwd("~/Box/Bernal_lab/Pinfish/ThermalStress_Project_data/Pinfish_ThermalStress_April2021/TagSeq_R_map_to_isoform/")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(DEGreport)
library(pheatmap)
library(RColorBrewer)
library(stringr)

allcountdata<-as.matrix(read.table("allcounts_isoform.txt"), sep="\t", row.names=1)
allcoldata<-read.table("pheno_data.txt")
allcoldata$treatment<-as.factor(allcoldata$treatment)

# Subset by tissue type
livercoldata<-subset(allcoldata, tissue == "liver")
livercountdata<-allcountdata[,c(1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49,51,53,55,57,59)]
musclecoldata<-subset(allcoldata, tissue == "muscle")
musclecountdata<-allcountdata[,c(2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60)]

# Check that things match
all(rownames(livercoldata) == colnames(livercountdata))
all(rownames(musclecoldata) == colnames(musclecountdata))

# Liver dataset
# Wald test for pairwise comparisons
Ldds<-DESeqDataSetFromMatrix(countData = livercountdata, colData = livercoldata, design =~ treatment)
Ldds$treatment <- relevel(Ldds$treatment, ref = "21")
Ldds<-DESeq(Ldds)

# Take a look at the results
Lres_27_v_24<-results(Ldds, alpha = 0.05, contrast = c("treatment", "27", "24"))
Lres_27_v_24
Lres_27_v_24_sigs<-subset(Lres_27_v_24, padj < 0.05)
summary(Lres_27_v_24)
# OUTPUT
#out of 122238 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 7, 0.0057%
#LFC < 0 (down)     : 4, 0.0033%
#outliers [1]       : 0, 0%
#low counts [2]     : 15, 0.012%
#(mean count < 0)

Lres_27_v_21<-results(Ldds, alpha = 0.05, contrast = c("treatment", "27", "21"))
Lres_27_v_21
Lres_27_v_21_sigs<-subset(Lres_27_v_21, padj < 0.05)
summary(Lres_27_v_21)
# OUTPUT
#out of 122238 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 693, 0.57%
#LFC < 0 (down)     : 1753, 1.4%
#outliers [1]       : 0, 0%
#low counts [2]     : 94499, 77%
#(mean count < 1)

Lres_24_v_21<-results(Ldds, alpha = 0.05, contrast = c("treatment", "24", "21"))
Lres_24_v_21
Lres_24_v_21_sigs<-subset(Lres_24_v_21, padj < 0.05)
summary(Lres_24_v_21)
#out of 122238 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 424, 0.35%
#LFC < 0 (down)     : 1311, 1.1%
#outliers [1]       : 0, 0%
#low counts [2]     : 101585, 83%
#(mean count < 2)

# LRT and PCA
# Plot just the significant genes on the PCA
liverLRT<-DESeq(Ldds, test="LRT", reduced=~1) # First run your likelihood ratio test
liverLRTres<-results(liverLRT, alpha = 0.05) # Generate your results object
summary(liverLRTres)
liverLRTvsd<-vst(liverLRT) # Transform the count data using variance stabilizing transformation
# We want just the significant genes, so isolate them here
liverLRTsigs<-subset(liverLRTres, padj < 0.05)
liver_sigs<-rownames(liverLRTsigs) # Isolate just the names of the significant genes
VSD.liver.subset <- liverLRTvsd[rownames(liverLRTvsd) %in% liver_sigs, ] # Extract transformed counts for significant genes
summary(VSD.liver.subset) # Make sure you've got them and they are correct
# Plot PCA of just significant genes, color by treatment
plotPCA(VSD.liver.subset, intgroup = "treatment")

# Muscle
Mdds<-DESeqDataSetFromMatrix(countData = musclecountdata, colData = musclecoldata, design =~ treatment)
Mdds$treatment <- relevel(Mdds$treatment, ref = "21")
Mdds<-DESeq(Mdds)

# Take a look at your results
Mres_27_v_24<-results(Mdds, alpha = 0.05, contrast = c("treatment", "27", "24"))
Mres_27_v_24
Mres_27_v_24_sigs<-subset(Mres_27_v_24, padj < 0.05)
summary(Mres_27_v_24)
# OUTPUT
#out of 100659 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 520, 0.52%
#LFC < 0 (down)     : 820, 0.81%
#outliers [1]       : 0, 0%
#low counts [2]     : 82635, 82%
#(mean count < 1)

Mres_27_v_21<-results(Mdds, alpha = 0.05, contrast = c("treatment", "27", "21"))
Mres_27_v_21
Mres_27_v_21_sigs<-subset(Mres_27_v_21, padj < 0.05)
summary(Mres_27_v_21)
# OUTPUT
#out of 100659 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 173, 0.17%
#LFC < 0 (down)     : 146, 0.15%
#outliers [1]       : 0, 0%
#low counts [2]     : 2, 0.002%
#(mean count < 0)

Mres_24_v_21<-results(Mdds, alpha = 0.05, contrast = c("treatment", "24", "21"))
Mres_24_v_21
Mres_24_v_21_sigs<-subset(Mres_24_v_21, padj < 0.05)
summary(Mres_24_v_21)
#out of 100659 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 941, 0.93%
#LFC < 0 (down)     : 435, 0.43%
#outliers [1]       : 0, 0%
#low counts [2]     : 82635, 82%
#(mean count < 1)

# LRT and PCA
# Plot just the significant genes on the PCA
musLRT<-DESeq(Mdds, test="LRT", reduced=~1) # First run your likelihood ratio test
musLRTres<-results(musLRT, alpha = 0.05) # Generate your results object
summary(musLRTres)
musLRTvsd<-vst(musLRT) # Transform the count data using variance stabilizing transformation
# We want just the significant genes, so isolate them here
musLRTsigs<-subset(musLRTres, padj < 0.05)
mus_sigs<-rownames(musLRTsigs) # Isolate just the names of the significant genes
VSD.mus.subset <- musLRTvsd[rownames(musLRTvsd) %in% mus_sigs, ] # Extract transformed counts for significant genes
summary(VSD.mus.subset) # Make sure you've got them and they are correct
# Plot PCA of just significant genes, color by treatment
plotPCA(VSD.mus.subset, intgroup = "treatment")

# Write output files to be used in GO term tests
M_27_v_21_df<-as.data.frame(Mres_27_v_21)
write.table(M_27_v_21_df, "muscle_27_21_pairwise_juv.tsv", sep = "\t", quote=F)
M_24_v_21_df<-as.data.frame(Mres_24_v_21)
write.table(M_24_v_21_df, "muscle_24_21_pairwise_juv.tsv", sep = "\t", quote=F)
M_27_v_24_df<-as.data.frame(Mres_27_v_24)
write.table(M_27_v_24_df, "muscle_27_24_pairwise_juv.tsv", sep = "\t", quote=F)

L_27_v_21_df<-as.data.frame(Lres_27_v_21)
write.table(L_27_v_21_df, "liver_27_21_pairwise_juv.tsv", sep = "\t", quote=F)
L_24_v_21_df<-as.data.frame(Lres_24_v_21)
write.table(L_24_v_21_df, "liver_24_21_pairwise_juv.tsv", sep = "\t", quote=F)
L_27_v_24_df<-as.data.frame(Lres_27_v_24)
write.table(L_27_v_24_df, "liver_27_24_pairwise_juv.tsv", sep = "\t", quote=F)
