##
## Title: KO filter for Heaney et al. 2021
## Author: Sanjeev V Namjoshi (snamjoshi87@utexas.edu)
## Date: February 23rd, 2017
## Revised: June 15th, 2019
##
# Load packages ---------------------------------
library(ggplot2)
library(reshape2)
library(biomaRt)
library(magrittr)
library(xlsx)
library(mgu74a.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(seqinr)
library(grid)
library(gridExtra)
library(edgeR)
library(limma)
library(Ckmeans.1d.dp)
library(pheatmap)
library(RColorBrewer)

# Functions ---------------------------------------------------------------
cutoff <- function(range) {
  
  above <- table(fullDist[range,"B_D_A_RipDist"]) +
    table(fullDist[range,"B_D_A_ParDist"]) +
    table(fullDist[range,"B_DDist"]) +
    table(fullDist[range,"D_A_ParDist"])
  
  
  below <- table(fullDist[length(range):nrow(ncAvg),"B_D_A_RipDist"]) +
    table(fullDist[length(range):nrow(ncAvg),"B_D_A_ParDist"]) +
    table(fullDist[length(range):nrow(ncAvg),"B_DDist"]) +
    table(fullDist[length(range):nrow(ncAvg),"D_A_ParDist"])
  
  return((above / (above + below))) 
}

anno <- function(x) {
  require(AnnotationDbi)
  require(org.Mm.eg.db)
  
  ### Annotate
  x$symbol <- mapIds(org.Mm.eg.db,
                     keys = row.names(x),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")
  
  ### Remove rows with no gene mapping
  x <- subset(x, select = c(17,1:16))
  x <- x[complete.cases(x), ]
  rownames(x) <- make.names(x$symbol, unique = TRUE)
  x[ ,1] <- NULL
  
  return(x)
}

# Load data ---------------------------------------------------------------
nc <- read.csv("finalCountData.csv", header = TRUE)
rownames(nc) <- nc[ ,1]
nc <- nc[ ,-1]

ncAvg <- data.frame(SAL_WT = rowMeans(nc[ ,1:3]),
                    SAL_KO = rowMeans(nc[ ,4:6]),
                    RO_WT = rowMeans(nc[ ,7:9]),
                    RO_KO = nc[ ,10],
                    RO_RAP_WT = rowMeans(nc[ ,11:13]),
                    RO_RAP_KO = rowMeans(nc[ ,14:16]))

ncAvg <- ncAvg[order(ncAvg$SAL_WT, decreasing = TRUE), ]

# Load other data  --------------------------------------------------------
### Brown data
brownDat <- read.xlsx("brown_list.xls", sheetIndex = 1, header = TRUE, startRow = 2)
brownProbes <- brownDat[ ,"Probe.Set"] %>% as.vector()
probes <- as.list(mgu74aALIAS2PROBE)
brownGenes <- probes[probes %in% brownProbes] %>% unlist() %>% names()

### Darnell data
darnellDat <- read.xlsx("darnell_list.xls", sheetIndex = 1, header = TRUE, startRow = 2)
darnellGenes <- darnellDat[ ,"Symbol.from.mm9"] %>% as.vector()

### Ascano PAR-CLIP data
ascanoParDat <- read.xlsx("ascano_par.xlsx", sheetIndex = 1, header = TRUE, startRow = 3, colIndex = 1:4)
ascanoParGenes <- ascanoParDat[ascanoParDat$Total.mRNA.binding.sites > 0, "Gene"] %>% as.vector()

### Ascano RIP-CHIP data
ascanoRipDat <- read.csv("ascano_rip.csv", header = TRUE, na.strings = "")
ascanoRipDat <- ascanoRipDat[ ,c(1,3)]
ascanoRipDat <- split(ascanoRipDat, ascanoRipDat$Target.bin)
ascanoRipDat <- ascanoRipDat$Target
ascanoRipGenes <- ascanoRipDat[ ,1] %>% as.vector()

### Miyashiro data
miyashiroDat <- read.csv("miyashiro.csv", header = TRUE)
miyashiroDat <- miyashiroDat[-c(1,1153:1159), ]
miyashiroDat <- miyashiroDat[ ,"UniGene"] %>% as.data.frame()
names(miyashiroDat) <- "unigene"
miyashiroDat$unigene <- as.character(miyashiroDat$unigene)
miyashiroDat$symbol <- mapIds(org.Hs.eg.db,
                              keys = miyashiroDat$unigene,
                              column = "SYMBOL",
                              keytype = "UNIGENE",
                              multiVals = "first")
miyashiroGenes <- miyashiroDat$symbol %>% na.omit() %>% as.vector()

### Consensus data
BDARIP <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 2, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1
BDAPAR <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 3, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1
BD <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 4, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1
DAPAR <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 5, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1

# Other data preload ------------------------------------------------------
### Consensus data
BDARIP <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 2, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1
BDAPAR <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 3, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1
BD <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 4, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1
DAPAR <- read.xlsx("consensus_FMRP.xlsx", sheetIndex = 5, header = FALSE, startRow = 2, stringsAsFactors = FALSE)$X1

brownGenes <- read.csv("brownGenes.csv", stringsAsFactors = FALSE)$x
darnellGenes <- read.csv("darnellGenes.csv", stringsAsFactors = FALSE)$x
ascanoParGenes <- read.csv("ascanoParGenes.csv", stringsAsFactors = FALSE)$x
ascanoRipGenes <- read.csv("ascanoRipGenes.csv", stringsAsFactors = FALSE)$x
miyashiroGenes <- read.csv("miyashiroGenes.csv", stringsAsFactors = FALSE)$x

# Visualize cutoff distributions ------------------------------------------
### Distributions across ranked saline genes for full data sets
fullDist <- data.frame(row.names = row.names(ncAvg) %>% toupper(),
                       order = 1:nrow(ncAvg),
                       brownDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (brownGenes %>% toupper()), "Brown et al. 2001", NA),
                       darnellDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (darnellGenes %>% toupper()), "Darnell et al. 2011", NA),
                       ascanoParDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (ascanoParGenes %>% toupper()), "Ascano et al. 2012 (PAR-CLIP)", NA),
                       ascanoRipDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (ascanoRipGenes %>% toupper()), "Ascano et al. 2012 (RIP-CHIP)", NA),
                       miyashiroRipDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (miyashiroGenes %>% toupper()), "Miyashiro et al. 2003", NA),
                       B_D_A_RipDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (BDARIP %>% toupper()), "B_D_A-RIP", NA),
                       B_D_A_ParDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (BDAPAR %>% toupper()), "B_D_A-PAR", NA),
                       B_DDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (BD %>% toupper()), "B_D", NA),
                       D_A_ParDist = ifelse((row.names(ncAvg) %>% toupper()) %in% (DAPAR %>% toupper()), "D_A-PAR", NA))

fullDist.m <- melt(fullDist, id.vars = "order")

fullDist.m$value <- factor(fullDist.m$value, levels = rev(c("Darnell et al. 2011", "Brown et al. 2001", "B_D_A-RIP", "B_D_A-PAR", "B_D", "D_A-PAR", "Ascano et al. 2012 (PAR-CLIP)", "Ascano et al. 2012 (RIP-CHIP)", "Miyashiro et al. 2003")))
fullDist.m$variable <- factor(fullDist.m$variable, levels = rev(c("darnellDist", "brownDist", "B_D_A_RipDist", "B_D_A_ParDist", "B_DDist", "D_A_ParDist", "ascanoParDist", "ascanoRipDist", "miyashiroRipDist")))

fullDist.m <- fullDist.m[complete.cases(fullDist.m), ]

# Ordered coverage matrix -------------------------------------------------
### Above/below cutoff line
coverage = NULL
for(i in 1:nrow(ncAvg)) {
  coverage[i] <- cutoff(1:i)
}

orderedCoverage <- data.frame(order = 1:nrow(ncAvg),
                              coverage = coverage * 100)

# Preload coverage --------------------------------------------------------
orderedCoverage <- read.csv("orderedCoverage.csv", stringsAsFactors = FALSE)

# Consensus cutoff --------------------------------------------------------
consensusCutoff <- orderedCoverage[orderedCoverage$coverage < 95, ] %>% nrow()

### Consensus percentage graph
a <- ggplot(orderedCoverage, aes(x = order, y = coverage)) +
  geom_line() +
  geom_vline(xintercept = consensusCutoff, color = "red", linetype = "dashed") +
  xlab("Ranked Saline Genes") +
  ylab("Percent Consensus \n Lists in Cutoff") +
  theme_bw()

a <- a + theme(axis.text.y = element_text(size = 12),
               axis.text.x = element_text(size = 12),
               text = element_text(size = 14))

### Consensus overlap graph
b <- ggplot(fullDist.m, aes(x = order, y = value)) +
  geom_point(alpha = 0.75, color = "dodgerblue") +
  geom_vline(xintercept = consensusCutoff, color = "red", linetype = "dashed") +
  xlab("") +
  ylab("") +
  theme_bw()

b <- b +   theme(axis.text.y = element_text(size = 12),
                 axis.text.x = element_text(size = 12))

### Fully arranged graphs (1234 x 546)
gA <- ggplotGrob(a)
gB <- ggplotGrob(b)
grid::grid.newpage()
grid::grid.draw(rbind(gB, gA))

# KO comparison to other data ---------------------------------------------
### Create enrichment ratios of WT over KO for each treatment
ncAvgDif <- data.frame(row.names = row.names(ncAvg),
                       SAL = ncAvg$SAL_WT / ncAvg$SAL_KO,
                       RO = ncAvg$RO_WT / ncAvg$RO_KO,
                       RO_RAP = ncAvg$RO_RAP_WT / ncAvg$RO_RAP_KO)

### Load targets and non-targets
targets <- c("Gabbr1", "Gabbr2", "Nisch", "Prrc2a", "Kif1a", "Grin2a", "Adgrb2", "Cyfip2", "Map1b", "Bsn", "Adcy1", "Pde2a", "Lingo1", "Grin2b")
nontargets <- c("Pabpc1", "Tmem65", "Hprt1", "St8sia3", "Sae1", "Glrb", "Gria2", "Tcerg1", "Eif3a", "Slc35f1", "Atp6ap2", "Vldlr")

### Create a data frame that includes the overlap of genes with our data set for each consensus gene list as well as the full Brown and Darnell lists
foldRange <- rbind(data.frame(group = "Brown et al. 2001", value = ncAvgDif[row.names(ncAvgDif) %in% brownGenes, ]$SAL),
                   data.frame(group = "Darnell et al. 2011", value = ncAvgDif[row.names(ncAvgDif) %in% darnellGenes, ]$SAL),
                   data.frame(group = "B_D_A-RIP", value = ncAvgDif[(row.names(ncAvgDif) %>% toupper()) %in% BDARIP, ]$SAL),
                   data.frame(group = "B_D_A-PAR", value = ncAvgDif[(row.names(ncAvgDif) %>% toupper()) %in% BDAPAR, ]$SAL),
                   data.frame(group = "B_D", value = ncAvgDif[(row.names(ncAvgDif) %>% toupper()) %in% BD, ]$SAL),
                   data.frame(group = "D_A-PAR", value = ncAvgDif[(row.names(ncAvgDif) %>% toupper()) %in% DAPAR, ]$SAL),
                   data.frame(group = "Non-Targets", value = ncAvgDif[row.names(ncAvgDif) %in% nontargets, ]$SAL))

### Determine average enrichment ratio range between WT and KO for all groups
foldRange.s <- split(foldRange, foldRange$group)

### Take the mean of the lower quartile of the ranges for all groups (indicated by the dotted line)
cutoff <- sapply(foldRange.s, function(x) quantile(x$value, prob = 0.25))[1:6] %>% mean()

# Univariate K-means for KO cutoff ----------------------------------------
ncAvgCut <- ncAvg[1:consensusCutoff, ]
ncAvgCut$"Saline" <- ncAvg[1:consensusCutoff, "SAL_WT"] / ncAvg[1:consensusCutoff, "SAL_KO"]
ncAvgCut$"Ro" <- ncAvg[1:consensusCutoff, "RO_WT"] / ncAvg[1:consensusCutoff, "RO_KO"]
ncAvgCut$"Ro+Rapa" <- ncAvg[1:consensusCutoff, "RO_RAP_WT"] / ncAvg[1:consensusCutoff, "RO_RAP_KO"]

### Univariate k-means clustering. Arrange output in data frame
clustDat <- data.frame(row.names = row.names(ncAvgCut),
                       enriched = ncAvgCut$Saline,
                       cluster = Ckmeans.1d.dp(ncAvgCut$Saline, 2)$cluster)

### Define a second cutoff by taking the max in the "background" cluster
cutoff2 <- split(clustDat, clustDat$cluster)[[1]]$enriched %>% max()
cutoff2

# Distributions above and below cutoff ------------------------------------
### Define knockout cutoff
KOcutoff <- mean(cutoff, cutoff2)
KOcutoff

### Filter data to obtain this cutoff
filteredTargets <- ncAvgCut[ncAvgCut$Saline > KOcutoff, 1:6]

### Create a filter variable that specifies signal or background
ncAvgCut$filter <- ifelse(ncAvgCut$Saline > KOcutoff, "Signal", "Background")

ncAvgCut[ncAvgCut$filter == "Signal", "Saline"] %>% mean() %>% round(3)
ncAvgCut[ncAvgCut$filter == "Signal", "Saline"] %>% sd() %>% round(3)
ncAvgCut[ncAvgCut$filter == "Background", "Saline"] %>% mean() %>% round(3)
ncAvgCut[ncAvgCut$filter == "Background", "Saline"] %>% sd() %>% round(3)


# Filter mitochonrial and glial cells -------------------------------------
mito <- read.csv("mito.csv", header = FALSE, stringsAsFactors = FALSE)$V1
glia <- read.csv("glial.csv", header = FALSE, stringsAsFactors = FALSE)$V1

ncAvgCutFilter <- ncAvgCut[!(row.names(ncAvgCut) %in% c(mito, glia)), ]
write.csv(ncAvgCutFilter, "counts_with_background.csv")

filteredTargets <- filteredTargets[!(row.names(filteredTargets) %in% c(mito, glia)), ]
finalCountData <- nc[row.names(nc) %in% row.names(filteredTargets), ]

# Filtration with limma ---------------------------------------------------
### Preprocess data columns
countData <- read.csv(file = "combinedCounts.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(countData) <- countData[ ,1]
countData[ ,1] <- NULL
cdMain <- subset(countData, select = c(7:12, 16:18, 3:6, 13:15)) %>% anno()
cdMain <- cdMain[ ,1:9]
cdMain <- cdMain[row.names(cdMain) %in% row.names(finalCountData), ]

### Normalize
dge <- DGEList(counts = cdMain)
dge <- calcNormFactors(dge)
dge <- cpm(dge, log = TRUE, prior.count = 3)

### Analyze with limma+voom
design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3)))
colnames(design) <- c("RO_RAP_WT", "RO_WT", "SAL_WT")

contrast.matrix <- makeContrasts(SAL_WT-RO_WT, RO_RAP_WT-RO_WT, levels = design)

v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
resultsRo <- topTable(fit2, coef = 1, number = 25000)
resultsRoRap <- topTable(fit2, coef = 2, number = 25000)

resultsRo$Threshold <- ifelse(resultsRo$adj.P.Val < 0.1, "Significant", "Non-Significant")
resultsRoRap$Threshold <- ifelse(resultsRo$adj.P.Val < 0.1, "Significant", "Non-Significant")









