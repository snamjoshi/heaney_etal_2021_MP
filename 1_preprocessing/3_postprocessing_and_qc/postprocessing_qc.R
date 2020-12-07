##
## Title: Postprocessing QC for Heaney et al. 2021
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)
## Date: February 23rd, 2017
## Revised: June 15th, 2019
##

# Load packages ---------------------------------
library(DESeq2)
library(magrittr)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)

# Load data ---------------------------------------------------------------
### Load count data
countData1 <- read.csv("countData_july.csv", header = TRUE)
rownames(countData1) <- countData1[ ,1]
countData1[ ,1] <- NULL

countData2 <- read.csv("countData_aug.csv", header = TRUE)
countData2[ ,1] <- NULL
names(countData2) <- gsub("_T1", "_T3", names(countData2))
names(countData2) <- gsub("_T2", "_T4", names(countData2))

countData <- cbind(countData1, countData2) %>% as.data.frame()

### Create metatable
metaTable <- data.frame(LibraryName = names(countData),
                        technical = c(rep(1:2, 21), rep(3:4, 21)))

metaTable$LibraryName <- as.factor(metaTable$LibraryName)
metaTable$technical <- as.factor(metaTable$technical)
metaTable$Counts <- paste(metaTable$LibraryName, "count", sep = ".")

rownames(metaTable) <- metaTable[ ,1]   # Set first column to row names
metaTable[ ,1] <- NULL

### Collapse technical replicates
dd <- DESeqDataSetFromMatrix(countData = countData,
                             colData = metaTable,
                             design = ~ technical)

dd$technical <- gsub("_T1|_T2|_T3|_T4", "", row.names(metaTable))
ddRep <- collapseReplicates(dd, groupby = dd$technical)
ddRep$technical <- as.factor(ddRep$technical)

### rlog transformation of data
rd <- rlog(ddRep)
ddRep <- estimateSizeFactors(ddRep)
rdCounts <- assay(rd) %>% as.data.frame()

# Count correlation -------------------------------------------------------
### Saline WT scatterplots
p2 <- ggplot(rdCounts, aes(SAL_WT1, SAL_WT2)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw() +

p3 <- ggplot(rdCounts, aes(SAL_WT2, SAL_WT3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### 546 x 441
# p3 + theme(axis.text.x = element_text(size = 12),
#            axis.text.y = element_text(size = 12),
#            text = element_text(size = 14))

p4 <- ggplot(rdCounts, aes(SAL_WT1, SAL_WT3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### Saline KO scatterplots
p5 <- ggplot(rdCounts, aes(SAL_KO1, SAL_KO2)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p6 <- ggplot(rdCounts, aes(SAL_KO2, SAL_KO3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p7 <- ggplot(rdCounts, aes(SAL_KO1, SAL_KO3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### Ro WT scatteplots
p8 <- ggplot(rdCounts, aes(RO_WT1, RO_WT2)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p9 <- ggplot(rdCounts, aes(RO_WT2, RO_WT3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p10 <- ggplot(rdCounts, aes(RO_WT1, RO_WT3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### Ro KO histogram
p11 <- ggplot(rdCounts, aes(RO_KO2, RO_KO2)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### Ro Rap WT scatterplots
p12 <- ggplot(rdCounts, aes(RO_RAP_WT1, RO_RAP_WT2)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p13 <- ggplot(rdCounts, aes(RO_RAP_WT2, RO_RAP_WT3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p14 <- ggplot(rdCounts, aes(RO_RAP_WT1, RO_RAP_WT3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### Ro Rap KO scatterplots
p15 <- ggplot(rdCounts, aes(RO_RAP_KO1, RO_RAP_KO2)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p16 <- ggplot(rdCounts, aes(RO_RAP_KO2, RO_RAP_KO3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

p17 <- ggplot(rdCounts, aes(RO_RAP_KO1, RO_RAP_KO3)) +
  geom_point(size = 0.5, alpha = 0.3) +
  theme_bw()

### Fully arranged (766 x 673)
grid.arrange(p2, p3, p4, p5,
             p6, p7, p8, p9,
             p10, p11, p12, p13,
             p14, p15, p16, p17,
             ncol = 4,
             nrow = 4)

# PCA ---------------------------------------------------------------------
# Arrange variables for plotting with ggplot2
rdNoTotal <- rd[ ,rd$technical %in% colData(rd)[3:18, 1]]
pcaDatNoTotal <- plotPCA(rdNoTotal, intgroup = "technical", returnData = TRUE)
percentVar <- round(100 * attr(pcaDatNoTotal, "percentVar"))

### 595 x 446
ggplot(pcaDatNoTotal, aes(PC1, PC2, color = technical, label = name)) +
  geom_point(size = 3) +
  geom_text(aes(label = name), size = 3, hjust = 0.3, vjust = -0.5) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  scale_color_manual(values = c(rep("darkolivegreen3", 1), rep("cadetblue3", 3), rep("dodgerblue3", 3), rep("chartreuse4", 3), rep("lightcoral", 3), rep("firebrick", 3))) +
  theme_bw() +
  xlab("Principal Component 1") +
  ylab("Principal Component 2") +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12)) +
  theme(legend.position = "none")

# Distance matrix ---------------------------------------------------------
### 632 x 446
sampleDists <- dist(t(assay(rdNoTotal)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(rdNoTotal)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)



