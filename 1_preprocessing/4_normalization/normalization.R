##
## Title: Normalization for Heaney et al. 2021
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)
## Date: February 23rd, 2017
## Revised: June 15th, 2019
##
# Load packages ---------------------------------
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(HTSFilter)
library(magrittr)
library(ggplot2)
library(reshape2)
library(ggthemes)

# Function ----------------------------------------------------------------
squash <- function(x) {
  output <- ((x - min(x)) / (max(x) - min(x))) * 100
  return(output)
}

percentChange <- function(new, old) {
  return(((new - old) / abs(old)) * 100)
}

# Count filtering ---------------------------------------------------------
countData <- read.csv(file = "combinedCounts.csv", header = TRUE)
rownames(countData) <- countData[ ,1]
countData[ ,1] <- NULL

countData$RO_KO1 <- countData$RO_KO2
countData$RO_KO3 <- countData$RO_KO2

cdMain <- subset(countData, select = c(16:18, 13:15, 10:12, 3, 22, 23, 7:9, 4:6))
conds <- sapply(names(cdMain), function(x) substr(x, 1, nchar(x)-1)) %>% unname

filter <- HTSFilter(x = cdMain, 
                    conds = conds,
                    normalization = "none")

cdFilter <- as.data.frame(filter$filteredData)
cdFilter[ ,c(11, 12)] <- NULL   # Remove duplicated Ro KO columns

# Annotation --------------------------------------------------------------
cdFilter$symbol <- mapIds(org.Mm.eg.db,
                          keys = row.names(cdFilter),
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

### Remove rows with no gene mapping
cdFilter <- subset(cdFilter, select = c(17, 1:16))
cdFilter <- cdFilter[complete.cases(cdFilter), ]
rownames(cdFilter) <- make.names(cdFilter$symbol, unique = TRUE)
cdFilter[ ,1] <- NULL

# Normalize counts --------------------------------------------------------

### Create metatable
metaTable <- data.frame(LibraryName = names(cdFilter),
                        LibraryLayout = rep("SINGLE", 16),
                        condition = c(rep("WT", 3), rep("KO", 3), rep("WT", 3), "KO", rep("WT", 3), rep("KO", 3)))
rownames(metaTable) <- metaTable[ ,1]   # Set first column to row names
metaTable[ ,1] <- NULL

### Get normalized counts
dd <- DESeqDataSetFromMatrix(countData = cdFilter,
                             colData = metaTable,
                             design = ~ 1)
dd <- estimateSizeFactors(dd)
nc <- as.data.frame(counts(dd, normalize = TRUE))

### Scale by BioAnalyzer values
BA <- read.csv("BA_values.csv", header = TRUE, stringsAsFactors = FALSE)
BA$mean <- rowMeans(BA[ ,2:4], na.rm = TRUE)

scale <- c(RO = BA$mean[2] / BA$mean[5],
           RO_RAP = BA$mean[3] / BA$mean[4],
           SAL = BA$mean[6] / BA$mean[7])

cdFilterScaled <- nc

cdFilterScaled[ ,10] <- cdFilterScaled[ ,10] / (scale[1] + 1)         # Scale Ro
cdFilterScaled[ ,14:16] <- cdFilterScaled[ ,14:16] / (scale[2] + 1)   # Scale Ro Rap
cdFilterScaled[ ,4:6] <- cdFilterScaled[ ,4:6] / (scale[3] + 1)       # Scale Sal

### Average across rows
nscMean <- data.frame(row.names = row.names(cdFilterScaled),
                      SAL_WT = rowMeans(cdFilterScaled[ ,1:3]),
                      SAL_KO = rowMeans(cdFilterScaled[ ,4:6]),
                      RO_WT = rowMeans(cdFilterScaled[ ,7:9]),
                      RO_KO = cdFilterScaled[ ,10],
                      RO_RAP_WT = rowMeans(cdFilterScaled[ ,11:13]),
                      RO_RAP_KO = rowMeans(cdFilterScaled[ ,14:16])) %>% squash()

# Visualize counts --------------------------------------------------------
### Create variables for plotting with ggplot2
nscMean.m <- melt(nscMean)
nscMean.m$Genotype <- c(rep("WT", 14543), rep("KO", 14543), rep("WT", 14543), rep("KO", 14543), rep("WT", 14543), rep("KO", 14543))
nscMean.m$variable <- c(rep("Saline", 29086), rep("Ro-25-6981", 29086), rep("Ro-25-6981 + Rapamycin", 29086))
nscMean.m$varOrder <- factor(nscMean.m$variable, levels = c("Saline", "Ro-25-6981", "Ro-25-6981 + Rapamycin"))

### Plot WT v. KO by treatment (1438 x 878)
ggplot(nscMean.m, aes(x = Genotype, y = value, color = Genotype)) +
  geom_point(position = position_jitter(width = 0.2, height = 0.2)) +
  facet_grid(~ varOrder) +
  xlab("") +
  ylab("Scaled Counts") +
  theme_few() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        text = element_text(size = 14))

# Visualize counts with FMRP targets --------------------------------------
### Load targets and non-targets
targets <- c("Gabbr1", "Gabbr2","Nisch", "Prrc2a", "Kif1a", "Grin2a", "Adgrb2", "Cyfip2", "Map1b", "Bsn", "Adcy1", "Pde2a", "Lingo1", "Grin2b")
nontargets <- c("Pabpc1", "Tmem65", "Hprt1", "St8sia3", "Sae1", "Glrb", "Gria2", "Tcerg1", "Eif3a", "Slc35f1", "Atp6ap2", "Vldlr")

### Prepare variables for plotting
nscRawTargets <- nscMean[row.names(nscMean) %in% c(targets, nontargets), ] %>% squash()
nscRawTargets$targets <- ifelse(row.names(nscRawTargets) %in% targets, "Target", "Non-target")
nscRawTargets.m <- melt(nscRawTargets)
nscRawTargets.m$Genotype <- c(rep("WT", 25), rep("KO", 25), rep("WT", 25), rep("KO", 25), rep("WT", 25), rep("KO", 25))
nscRawTargets.m$variable <- c(rep("Saline", 50), rep("Ro-25-6981", 50), rep("Ro-25-6981 + Rapamycin", 50))
nscRawTargets.m$group <- rep(1:25, 6)
nscRawTargets.m$varOrder <- factor(nscRawTargets.m$variable, levels = c("Saline", "Ro-25-6981", "Ro-25-6981 + Rapamycin"))
nscRawTargets.m$GenotypeOrder <- factor(nscRawTargets.m$Genotype, levels = c("WT", "KO"))

### Plot WT v. KO for FMRP targets and non-targets (1438 x 878)
ggplot(nscRawTargets.m, aes(x = GenotypeOrder, y = value, group = group, color = targets)) +
  geom_point() +
  geom_line() +
  facet_grid(~ varOrder) +
  xlab("") +
  ylab("Scaled Counts") +
  theme_few() +
  scale_color_discrete(name = "FMRP Targets") +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        text = element_text(size = 14))

# Scaled counts as percent change -----------------------------------------
new <- nscRawTargets[nscRawTargets$targets == "Target", 1:6] %>% colMeans()
old <- nscRawTargets[nscRawTargets$targets == "Non-target", 1:6] %>% colMeans()

percentChange(new[1], old[1]) %>% round(2)
percentChange(new[3], old[3]) %>% round(2)
percentChange(new[5], old[5]) %>% round(2)






