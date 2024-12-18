setwd("~/Box/Bernal_lab/Pinfish/ThermalStress_Project_data/Pinfish_ThermalStress_December_2020/TagSeq_R_map_to_isoform/")

library(ggplot2)
library(DESeq2)
library(dplyr)
library(DEGreport)
library(pheatmap)
library(RColorBrewer)
library(stringr)

adultallcountdata<-as.matrix(read.table("allcounts_isoform.txt"), sep="\t", row.names=1)
# Remove fish that should be excluded based on size
adultallcountdata<-subset(adultallcountdata, select=-c(C1.3L.trim.sam.counts.isoform, C1.3M.trim.sam.counts.isoform, C2.4L.trim.sam.counts.isoform, C2.4M.trim.sam.counts.isoform, H1.4L.trim.sam.counts.isoform, H1.4M.trim.sam.counts.isoform, H2.4L.trim.sam.counts.isoform, H2.4M.trim.sam.counts.isoform, M1.2L.trim.sam.counts.isoform, M1.2M.trim.sam.counts.isoform))
adultallcoldata<-read.table("phenotype_data.txt")
adultallcoldata$treatment<-as.factor(adultallcoldata$treatment)
adultallcoldata <- adultallcoldata %>%
  mutate(weight_bin = if_else(
    weight > 20, "large", "small"
  ))
adultallcoldata<-subset(adultallcoldata, weight_bin == "large")
adultallcoldata$exposure_bin<-c("short","short","long","long","long","long",
                                "long",'long','long','long','long','long',
                                'short','short','short','short','short',
                                'short','short','short','short','short',
                                'short','short',"long",'long','long','long','long','long',
                                'medium','medium','medium','medium')

# Check columns and rows have the same name
all(rownames(adultallcoldata) == colnames(adultallcountdata))

# Subset by tissue type
musclecoldata<-subset(adultallcoldata, tissue == "muscle")
musclecountdata<-subset(adultallcountdata, select=c("C1.1M.trim.sam.counts.isoform","C1.5M.trim.sam.counts.isoform","C1.6M.trim.sam.counts.isoform","C2.2M.trim.sam.counts.isoform","C2.5M.trim.sam.counts.isoform","C2.7M.trim.sam.counts.isoform","H1.1M.trim.sam.counts.isoform","H1.2M.trim.sam.counts.isoform","H1.3M.trim.sam.counts.isoform","H2.1M.trim.sam.counts.isoform","H2.2M.trim.sam.counts.isoform","H2.3M.trim.sam.counts.isoform","M1.1M.trim.sam.counts.isoform","M1.3M.trim.sam.counts.isoform","M1.4M.trim.sam.counts.isoform","M2.2M.trim.sam.counts.isoform","M2.4M.trim.sam.counts.isoform"))
all(rownames(musclecoldata) == colnames(musclecountdata))

livercoldata<-subset(adultallcoldata, tissue == "liver")
livercountdata<-subset(adultallcountdata, select=c("C1.1L.trim.sam.counts.isoform","C1.5L.trim.sam.counts.isoform","C1.6L.trim.sam.counts.isoform","C2.2L.trim.sam.counts.isoform","C2.5L.trim.sam.counts.isoform","C2.7L.trim.sam.counts.isoform","H1.1L.trim.sam.counts.isoform","H1.2L.trim.sam.counts.isoform","H1.3L.trim.sam.counts.isoform","H2.1L.trim.sam.counts.isoform","H2.2L.trim.sam.counts.isoform","H2.3L.trim.sam.counts.isoform","M1.1L.trim.sam.counts.isoform","M1.3L.trim.sam.counts.isoform","M1.4L.trim.sam.counts.isoform","M2.2L.trim.sam.counts.isoform","M2.4L.trim.sam.counts.isoform"))
all(rownames(livercoldata) == colnames(livercountdata))

# Run a likelihood ratio test to determine if the effect of exposure time is significant
# First muscle
Mdds<-DESeqDataSetFromMatrix(countData = musclecountdata, colData = musclecoldata, design =~ resp + exposure_bin + treatment)
Mdds<-DESeq(Mdds)

# This gives you the number of genes due to exposure time
musLRTexp<-DESeq(Mdds, test="LRT", reduced=~resp + treatment)
musLRTexp.res<-results(musLRTexp, alpha = 0.05)
summary(musLRTexp.res)
# Zero genes significant due to exposure time, okay to exclude

# Then liver
Ldds<-DESeqDataSetFromMatrix(countData = livercountdata, colData = livercoldata, design =~ resp + exposure_bin + treatment)
Ldds<-DESeq(Ldds)

# This gives you the number of genes significant due to exposure time
liverLRTexp<-DESeq(Ldds, test = "LRT", reduced =~ resp + treatment)
liverLRTexp.res<-results(liverLRTexp, alpha = 0.05)
summary(liverLRTexp.res)
# Two genes significant due to exposure time, okay to exclude

# The code from above shows that the effect of exposure time is negligible and it is okay to exclude it from our further analyses

# Now we can do pairwise comparisons for liver and muscle

## MUSCLE ##
Mdds<-DESeqDataSetFromMatrix(countData = musclecountdata, colData = musclecoldata, design =~ resp + treatment)
Mdds<-DESeq(Mdds)

Mres_27_v_24 <- results(Mdds, alpha = 0.05, contrast = c("treatment","27","24"))
Mres_27_v_24
Mres_27_v_24_sigs<- subset(Mres_27_v_24, padj < 0.05)
Mres_27_v_24_sigs
summary(Mres_27_v_24)
# OUTPUT:
#out of 83253 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2, 0.0024%
#LFC < 0 (down)     : 6, 0.0072%
#outliers [1]       : 100, 0.12%
#low counts [2]     : 3069, 3.7%
#(mean count < 0)

Mres_27_v_21 <- results(Mdds, alpha = 0.05, contrast = c("treatment", "27", "21"))
Mres_27_v_21
Mres_27_v_21_sigs <- subset(Mres_27_v_21, padj < 0.05)
summary(Mres_27_v_21)
# OUTPUT:
#out of 83253 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 41, 0.049%
#LFC < 0 (down)     : 48, 0.058%
#outliers [1]       : 100, 0.12%
#low counts [2]     : 76896, 92%
#(mean count < 12)

Mres_24_v_21 <- results(Mdds, alpha = 0.05, contrast=c("treatment", "24", "21"))
Mres_24_v_21
Mres_24_v_21_sigs<- subset(Mres_24_v_21, padj < 0.05)
Mres_24_v_21_sigs
summary(Mres_24_v_21)
# OUTPUT:
#out of 83253 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 7, 0.0084%
#LFC < 0 (down)     : 6, 0.0072%
#outliers [1]       : 100, 0.12%
#low counts [2]     : 50258, 60%
#(mean count < 0)

# LRT and PCA
musLRT<-DESeq(Mdds, test="LRT", reduced=~resp)
musLRTres<-results(musLRT, alpha = 0.05)
summary(musLRTres)
musLRTvsd<-vst(musLRT)
musLRTsigs<-subset(musLRTres, padj < 0.05)
mus_sigs<-rownames(musLRTsigs)
VSD.mus.subset <- musLRTvsd[rownames(musLRTvsd) %in% mus_sigs, ]
summary(VSD.mus.subset)
plotPCA(VSD.mus.subset, intgroup = "treatment")
adult_muscle_PCA_dim<-plotPCA(VSD.mus.subset, intgroup = "treatment", returnData = T)
ggplot(adult_muscle_PCA_dim, aes(PC1, PC2)) + geom_point(color = as.factor(adult_muscle_PCA_dim$treatment)) +
  labs(x = "PC1: 31% variance", y = "PC2: 23% variance") + ggtitle("Adult muscle")
ggsave("~/Box/Bernal_lab/Pinfish/ThermalStress_Project_PAPER/Figures/Figure4_adult_muscle.svg", width = 5, height = 3)

# Now run the same for liver
# LIVER
Ldds<-DESeqDataSetFromMatrix(countData = livercountdata, colData = livercoldata, design =~ resp + treatment)
Ldds$treatment<-relevel(Ldds$treatment, ref="21")
Ldds<-DESeq(Ldds)

Lres_27_v_24 <- results(Ldds, alpha = 0.05, contrast = c("treatment","27","24"))
Lres_27_v_24
Lres_27_v_24_sigs<- subset(Lres_27_v_24, padj < 0.05)
Lres_27_v_24_sigs
summary(Lres_27_v_24)
# OUTPUT:
#out of 115412 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 10, 0.0087%
#LFC < 0 (down)     : 19, 0.016%
#outliers [1]       : 133, 0.12%
#low counts [2]     : 3760, 3.3%
#(mean count < 0)

Lres_27_v_21 <- results(Ldds, alpha = 0.05, contrast = c("treatment","27","21"))
Lres_27_v_21
Lres_27_v_21_sigs<- subset(Lres_27_v_21, padj < 0.05)
Lres_27_v_21_sigs
summary(Lres_27_v_21)
# OUTPUT:
#out of 115412 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 19, 0.016%
#LFC < 0 (down)     : 24, 0.021%
#outliers [1]       : 133, 0.12%
#low counts [2]     : 24480, 21%
#(mean count < 0)

Lres_24_v_21 <- results(Ldds, alpha = 0.05, contrast = c("treatment","24","21"))
Lres_24_v_21
Lres_24_v_21_sigs<- subset(Lres_24_v_21, padj < 0.05)
Lres_24_v_21_sigs
summary(Lres_24_v_21)
# OUTPUT:
#out of 115412 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 26, 0.023%
#LFC < 0 (down)     : 12, 0.01%
#outliers [1]       : 133, 0.12%
#low counts [2]     : 19837, 17%
#(mean count < 0)

# LRT and PCA
livLRT<-DESeq(Ldds, test="LRT", reduced=~resp)
livLRTres<-results(livLRT, alpha = 0.05)
summary(livLRTres)
livLRTvsd<-vst(livLRT)
# We want just the significant genes
livLRTsigs<-subset(livLRTres, padj < 0.05)
liv_sigs<-rownames(livLRTsigs)
VSD.liv.subset <- livLRTvsd[rownames(livLRTvsd) %in% liv_sigs, ]
summary(VSD.liv.subset)
plotPCA(VSD.liv.subset, intgroup = "treatment") # Here a bit of overlap between control and +3, but clear distinction between both of those and +6
adult_liver_PCA_dim<-plotPCA(VSD.liv.subset, intgroup = "treatment", returnData = T)
ggplot(adult_liver_PCA_dim, aes(PC1, PC2)) + geom_point(color = as.factor(adult_liver_PCA_dim$treatment)) +
  labs(x = "PC1: 38% variance", y = "PC2: 27% variance") + ggtitle("Adult liver")
ggsave("~/Box/Bernal_lab/Pinfish/ThermalStress_Project_PAPER/Figures/Figure4_adult_liver.svg", width = 5, height = 3)

# Write output files to be used in GO term tests
M_27_v_21_df<-as.data.frame(Mres_27_v_21)
write.table(M_27_v_21_df, "muscle_27_21_pairwise_adult.tsv", sep = "\t", quote=F)
M_24_v_21_df<-as.data.frame(Mres_24_v_21)
write.table(M_24_v_21_df, "muscle_24_21_pairwise_adult.tsv", sep = "\t", quote=F)
M_27_v_24_df<-as.data.frame(Mres_27_v_24)
write.table(M_27_v_24_df, "muscle_27_24_pairwise_adult.tsv", sep = "\t", quote=F)

L_27_v_21_df<-as.data.frame(Lres_27_v_21)
write.table(L_27_v_21_df, "liver_27_21_pairwise_adult.tsv", sep = "\t", quote=F)
L_24_v_21_df<-as.data.frame(Lres_24_v_21)
write.table(L_24_v_21_df, "liver_24_21_pairwise_adult.tsv", sep = "\t", quote=F)
L_27_v_24_df<-as.data.frame(Lres_27_v_24)
write.table(L_27_v_24_df, "liver_27_24_pairwise_adult.tsv", sep = "\t", quote=F)