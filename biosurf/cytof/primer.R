## Installing the R packages you will need

## From CRAN
install.packages(c("scales", "RColorBrewer", "MASS", "Rtsne", "kohonen", "miscTools", "gplots", "Radviz", "igraph"))

## From BioConductor
source("https://bioconductor.org/biocLite.R")
biocLite(c("flowCore", "cytofkit", "ConsensusClusterPlus", "cydar", "flowCL", "CATALYST", "ncdfFlow", "edgeR"))

## From source
# cytofCore
library(devtools)
install_github("nolanlab/cytofCore")

# MEM
library(devtools)
download.file(url = "http://www.nature.com/nmeth/journal/v14/n3/extref/nmeth.4149-S5.zip", destfile = "./MEM.zip")
unzip("MEM.zip")
setwd("./MEM")
build()
install.packages("./", repos = NULL, type="source")
unlink("../MEM/", recursive = TRUE) # optional cleanup


## The Flow Cytometry Standard (FCS) format

# Load the flowCore package
library(flowCore)

# Set working directory (directory where you placed the FCS files)
setwd("~/Dropbox/Research/Projects/Ongoing/Stanford/CyTOF/data/HIMC_controls/raw/")
# setwd("~/cytof_primer/fcs_files/")

# Read an fcs file
fcs_files <- list.files(pattern='.fcs$', full=TRUE, ignore.case = TRUE)
fcs <- read.FCS(filename=fcs_files[1], transformation=FALSE)

# Extract expression matrix
exprs <- fcs@exprs

# Explore the channels (columns in the expression matrix)
fcs@parameters@data

# Make colnames human readable using information in the parameter data slot
markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
colnames(exprs)[which(!is.na(markers))] <- markers[which(!is.na(markers))]

# Merging fcs files (in case of interrupted runs)
library(cytofCore)
dir.create("./combined")
cytofCore.concatenateDirectoryFiles(inputDir="./",outputDir="./combined",pattern=NULL,overwrite=F,timeParam="time")

# CyTOF data preprocessing
pregating_channels <- c("Bead", "DNA1", "DNA2", "Dead")
lineage_channels <- c("CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25")
instrument_channels <- c("Event_length", "Time", "Event_length", "Center", "Offset", "Width", "Residual")

## Data transformations

### Randomization

# Check if counts are randomized (expression matrix should only contain integers)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
table(is.wholenumber(exprs[,c(pregating_channels, lineage_channels)]))

# You may revert back to original counts (only applicable to measured channels - i.e. exclude time, event length, etc.)
exprs <- cbind(ceiling(exprs[,c(pregating_channels, lineage_channels)]), exprs[,instrument_channels])

### ArcSinh transformations
library(MASS)
library(RColorBrewer)
k <- 11; my.cols <- rev(brewer.pal(k, "RdYlBu"))

# Set co-factor
cofac <- 2

# Make arcsinh transformed expression matrix (with the exception of time and event_length, which should remain linear)
exprs_trans <- cbind(asinh(exprs[,c(pregating_channels, lineage_channels)]/cofac), exprs[,instrument_channels])

# Plot two-marker example of transformation effect
plot(exprs[,c("CD8", "CD4")], pch=".", col="grey", main="CD8 vs CD4")
z <- kde2d(exprs[,"CD8"], exprs[,"CD4"], n=50)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

plot(exprs_trans[,c("CD8", "CD4")], pch=".", col="grey", main="CD8 vs CD4 (transformed counts)")
z <- kde2d(exprs_trans[,"CD8"], exprs_trans[,"CD4"], n=50)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)


## Bead Normalization

# CODE #

## Data transformations

## Pre-gating

# For the next step, set working directory to the normalized samples
setwd("~/Dropbox/Research/Projects/Ongoing/Stanford/CyTOF/data/HIMC_controls/normed/")
# setwd("~/cytof_primer/fcs_files/normed/")

# Read in your normalized files and concatenate in a data frame (in this example we subset each matrix to 10,000 random events)
fcs_files <- list.files(pattern='.fcs$', full=TRUE, ignore.case = TRUE)
exprs_set <- data.frame()
sample <- c()
for(i in fcs_files) {
  fcs <- read.FCS(filename=i, transformation=FALSE) 
  exprs <- fcs@exprs
  markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
  colnames(exprs)[which(!is.na(markers))] <- markers[which(!is.na(markers))]
  exprs_set <- rbind(exprs_set, exprs[sample(nrow(exprs), 10000),])
  sample <- append(sample, rep(i, 10000))
}
exprs_set$sample <- sample
instrument_channels <- append(instrument_channels, "sample")

# Make arcsinh transformed expression matrix (excluding channels that should remain linear - e.g. time, event length, etc.)
cofac <- 2 # For v1 and v2 use a co-factor of 5, for Helios use a co-factor of 2
exprs_set_trans <-cbind(asinh(exprs_set[,c(pregating_channels, lineage_channels)]/cofac), exprs_set[,instrument_channels])

## Batch correction

# visualize densities of pre-gating channels in each sample
library(scales)
library(RColorBrewer)

cols <- brewer.pal(5,'Set1')
par(mfrow=c(1,5))

for(c in 1:length(pregating_channels)) {
  plot(1, type = "n", xlim = c(0, max(exprs_set_trans[,pregating_channels[c]])), ylim = c(1, length(fcs_files) + 2), axes = FALSE, xlab = "", ylab = "", main=pregating_channels[c], cex.main=3)
  for(s in 1:length(fcs_files)){
    dens <- density(exprs_set_trans[exprs_set_trans$sample==fcs_files[s],pregating_channels[c]])
    dens$y = dens$y+(length(fcs_files)-s+1)
    par(bg=NA)
    lines(dens)
    polygon(dens, col=alpha(cols[c],0.4), border = NA)
    abline(h =(length(fcs_files)-s+1), lwd = 0.5)
  }
}

dev.off()

# Create two toy batches of bimodal densities (the cydar packages only takes lists of lists as input)
batches <- list()
x <- list()
ex <- matrix(c(rnorm(500,2,1), rnorm(500,25,2)), ncol=1)
colnames(ex) <- "CDx"
x[["sample_Y"]] <- ex
batches[["Batch1"]] <- x
ex <- matrix(c(rnorm(500,2.1,1.1), rnorm(500,15,2)), ncol=1)
colnames(ex) <- "CDx"
x[["sample_Y"]] <- ex
batches[["Batch2"]] <- x

# Plot the marker distributions
plot(density(batches$Batch1$sample_Y), main="uncorrected samples", ylim=c(0,0.09))
points(density(batches$Batch2$sample_Y), type="l", col="red")

# Run the normalization
library(cydar)
corrected <- normalizeBatch(batches, batch.comp = NULL, mode="warp")

# Plot the results
plot(density(corrected$Batch1$sample_Y), main="batch corrected samples")
points(density(corrected$Batch2$sample_Y), type="l", col="red")


### Gating for cells (beads vs. DNA)

# Load necessary packages
library(MASS)
library(RColorBrewer)

# Make color palette using color brewer
k <- 11; my.cols <- rev(brewer.pal(k, "RdYlBu"))

sub <- exprs_set_trans[exprs_set_trans$sample==fcs_files[4],c("Bead", "DNA1")]
z <- kde2d(x=sub[,1], y=sub[,2], n=50, h=max(bandwidth.nrd(sub[,1]), bandwidth.nrd(sub[,2])))

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="grey", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 0; right <- 3; lower <- 6; upper <- 9
rect(left,lower,right,upper,border="red", lwd=3)
cells <- exprs_set_trans[exprs_set_trans[,"Bead"] < right,]; cells <- cells[cells[,"DNA1"] > lower & cells[,"DNA1"] < upper,]
paste("cells: ", round((nrow(cells)/nrow(exprs_set_trans))*100, digits=2), "%", sep="")

# Gate for intact cells (DNA1 vs DNA2)
sub <- cells[,c("DNA1", "DNA2")]
z <- kde2d(sub[,1], sub[,2], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2])
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 6.4; right <- 8.1; lower <- 7; upper <- 8.7
rect(left,lower,right,upper,border="red", lwd=3)
intact <- cells[cells[,"DNA1"] > left & cells[,"DNA1"] < right,]; intact <- intact[intact[,"DNA2"] > lower & intact[,"DNA2"] < upper,]
paste("intact cells: ", round((nrow(intact)/nrow(cells))*100, digits=2), "%", sep="")

# Gate for singlets (event length vs. DNA)
sub <- intact[,c("Event_length", "DNA1")]
z <- kde2d(sub[,1], sub[,2], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 11; bottom <- min(sub[,2]); right <- 22; upper <- max(sub[,2])
rect(left,bottom,right,upper,border="red", lwd=3) # Draw cells gate (format: ll x,y; ur x,y)
singlets <- intact[intact[,"Event_length"] > left & intact[,"Event_length"] < right,]
paste("Intact singlets: ", round((nrow(singlets)/nrow(intact))*100, digits=2), "%", sep="")

# Gate for live
sub <- singlets[,c("Dead", "DNA1")]
z <- kde2d(sub[,1], sub[,2], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 0; bottom <- min(sub[,2]); right <- 5; upper <- max(sub[,2])
rect(left,bottom,right,upper,border="red", lwd=3)
live <- singlets[singlets[,"Dead"] < right,]
paste("Live intact singlets: ", round((nrow(live)/nrow(singlets))*100, digits=2), "%", sep="")

for (i in fcs_files) {
  print(paste("Total recovery for ", i, ": ", round((nrow(live[live$sample==i,])/nrow(exprs_set_trans[exprs_set_trans$sample==i,]))*100, digits=2), "%", sep=""))
}

# Save FCS containing live intact singlets
ff2 <- flowFrame(xt_live_intact_singlets, ff@parameters, ff@description)
suppressWarnings(write.FCS(ff2, filename = paste("./",  file, "live_intact_singlets", sep=""), what="numeric", delimiter = "\\"))


#### LOAD FILE
# save(live, lineage_channels, instrument_channels, pregating_channels, fcs_files, file = "~/Dropbox/Research/Papers/Ongoing/CyTOF primer/CyTOF_wiki/live.Rdata")
load("~/Dropbox/Research/Papers/Ongoing/CyTOF primer/CyTOF_wiki/live.Rdata")


## Dimensionality reduction

### principal component analysis
pca <- prcomp(live[live$sample==fcs_files[1],lineage_channels], scale. = TRUE)
plot(pca$x[,1:2], pch=21, bg=alpha("slategray4",0.6), col="grey")

### t-SNE
# Remove duplicates
table(duplicated(live[live$sample==fcs_files[1],lineage_channels]))
# No duplicates here, but if so, they can be removed:
live <- live[-which(duplicated(live[live$sample==fcs_files[1],lineage_channels]))]

library(Rtsne)
set.seed(42)
tsne <- Rtsne(live[live$sample==fcs_files[1],lineage_channels])
plot(tsne$Y, pch=21, bg=alpha("slategray4",0.6), col="grey", xlab="", ylab="", lwd=0.5)

# Color by marker expression
par(mfrow=c(2,2))
for(i in c("CD4", "CD8", "CD19", "CD56")) {
  cols_ramp <- colorRamp(c("slategray4", "red"))
  marker_exprs_scaled <- (live[live$sample==fcs_files[1],i] - min(live[live$sample==fcs_files[1],i]))/diff(range(live[live$sample==fcs_files[1],i]))
  cols <- rgb(cols_ramp(marker_exprs_scaled), maxColorValue = 256)
  plot(tsne$Y, pch=21, bg=alpha(cols,0.6), col="grey", xlab="", ylab="", lwd=0.5, main=i)
}

## Cell subset detection (clustering)

### Phenograph
library(cytofkit)
clusters_pg <- cytof_cluster(xdata = live[live$sample==fcs_files[1],lineage_channels], method = "Rphenograph")

# Visualize results on PCA
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
plot(pca$x[,1:2], pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey")

# Visualize results on t-SNE plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
plot(tsne$Y, pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey", xlab="", ylab="", lwd=0.5)

# Sample 5000 random events from each sample
sub <- NULL
for(i in 1:4) {
  temp <- live[live$sample==fcs_files[i],]
  sub <- rbind(sub, temp[sample(nrow(temp), 5000, replace=FALSE),])
}
# Check for duplicated events (remove if applicable)
table(duplicated(sub[,lineage_channels]))
sub <- sub[-which(duplicated(sub[,lineage_channels]))]
# Cluster using phenograph
clusters_pg <- cytof_cluster(xdata = sub[,lineage_channels], method = "Rphenograph")
# Run t-SNE
set.seed(42)
tsne <- Rtsne(sub[,lineage_channels])
# Plot each sample on t-SNE canvas consisting of all events
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
par(mfrow=c(2,2))
for(i in 1:4) {
  plot(tsne$Y, pch=19, col="slategray4", xlab="", ylab="", lwd=0.5, main=paste("sample", i))
  points(tsne$Y[sub$sample==unique(sub$sample)[i],1:2], pch=21, bg=alpha(col=cols[clusters_pg[sub$sample==unique(sub$sample)[i]]],0.6), col="grey", xlab="", ylab="", lwd=0.5)
}

### FlowSOM clustering
library(cytofkit)
clusters_fs <- cytof_cluster(xdata = live[live$sample==fcs_files[1],lineage_channels], method = "FlowSOM", FlowSOM_k=25)

# Visualize results on t-SNE plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_fs))
plot(tsne$Y, pch=21, bg=alpha(col=cols[clusters_fs],0.6), col="grey", xlab="", ylab="", lwd=0.5)

### SPADE clustering

## Visualizaing expression in different clusters

### Simple heatmap
library(miscTools)
sub_matrix <- live[live$sample==fcs_files[1],lineage_channels]
cluster_matrix <- NULL
for(i in 1:max(clusters_pg)) {
  cluster_matrix <- rbind(cluster_matrix, colMedians(sub_matrix[clusters_pg==i,]))
}
library(gplots)
cols = colorRampPalette(c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[1]))
par(mar = c(2,2,2,2))
heatmap.2(t(cluster_matrix), col=cols, trace="none", density.info = "none", sepcolor = "white", sepwidth = c(0.001, 0.001), colsep=c(1:ncol(t(cluster_matrix))), rowsep=c(1:nrow(t(cluster_matrix))), xlab="cluster", ylab="channel", scale="row")

cluster_matrix_pseudoMEM <- NULL
for(i in 1:max(clusters_pg)) {
  cluster_matrix_pseudoMEM <- rbind(cluster_matrix_pseudoMEM, colMedians(sub_matrix[clusters_pg==i,]) - median(as.vector(as.matrix(sub_matrix[clusters_pg!=i,]))))
}
heatmap.2(t(cluster_matrix_pseudoMEM), col=cols, trace="none", density.info = "none", sepcolor = "white", sepwidth = c(0.001, 0.001), colsep=c(1:ncol(t(cluster_matrix))), rowsep=c(1:nrow(t(cluster_matrix))), xlab="phenograph cluster", ylab="channel")

### MEM

library(MEM)
# Calculate enrichment scores
mem <- MEM(cbind(live[live$sample==fcs_files[1],lineage_channels], cluster<-clusters_pg))
mem <- round(mem[[5]][[1]])
mem <- mem[order(as.numeric(rownames(mem))),]
threshold <- 4
labels <- list()
for(i in 1:nrow(mem)) {
  lab <- mem[i,abs(mem[i,])>threshold]
  lab[lab>0] <- paste("+",lab[lab>0], sep="")
  lab <- paste(names(lab), lab, sep="")
  labels[[rownames(mem)[i]]] <- lab
}

# MEM comes with a heatmap function based on heatmap.2, but if you wish to customize appearance, use the following

# Generating "human readable" labels for the clusters as discussed in the Nat. Meth. article require some extra coding
# First set threshold enrichment score for including a marker in label (CDx+/-cutoff)
cutoff = 2
mem_labels <- round(mem$MEM_matrix[[1]])
paste(names(mem_labels[1,][abs(mem_labels[1,])>0]), sort(mem_labels[1,][abs(mem_labels[1,])>0],decreasing = TRUE), sep="")
# Plot heatmap
cols = colorRampPalette(c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[9], brewer.pal(9, "Set1")[1]))
par(mar = c(2,2,2,2))
heatmap.2(t(mem$MEM_matrix[[1]]), col=cols, trace="none", density.info = "none", sepcolor = "white", sepwidth = c(0.001, 0.001), colsep=c(1:ncol(t(mem$MEM_matrix[[1]]))), rowsep=c(1:nrow(t(mem$MEM_matrix[[1]]))), xlab="phenograph cluster", ylab="channel")


### Radviz
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
library(Radviz)
cell.S <- make.S(colnames(sub_matrix))
cell.sim <- cosine(as.matrix(sub_matrix))
optim.cell <- do.optim(cell.S,cell.sim,iter=100,n=1000)
cell.S <- make.S(tail(optim.cell$best,1)[[1]])
cell.rv <- do.radviz(as.matrix(sub_matrix),cell.S)
layout(matrix(c(1,1,1,1,1,2), nrow=1,ncol=6))
plot(cell.rv, point.shape=19, point.color=alpha(cols[clusters_pg], 0.6))
plot.new()
legend("center", legend=paste("cluster", 1:max(clusters_pg)), col=cols, cex=1, pch=16, bty='n')

pop.rv <- do.radviz(cluster_matrix,cell.S)
size <- table(clusters_pg)
layout(matrix(c(1,1,1,1,1,2), nrow=1,ncol=6))
bubbleRadviz(pop.rv, bubble.color=alpha(cols[1:nrow(cluster_matrix)],0.8), bubble.size=log(size[1:nrow(cluster_matrix)]), scale=0.2, decreasing=TRUE)
plot.new()
legend("center", legend=paste("cluster", 1:max(clusters_pg)), col=cols, cex=1, pch=16, bty='n')

## Differential abundance of cell populations and proteins

# Import data
library(flowCore)
setwd("~/Dropbox/Research/Papers/Ongoing/CyTOF primer/CyTOF_wiki/cancer_data/")
fcs_files <- list.files(pattern='.fcs$', full=TRUE, ignore.case = TRUE)
exprs_set <- data.frame()
sample <- c()
for(i in fcs_files) {
  fcs <- read.FCS(filename=i, transformation=FALSE) 
  exprs <- fcs@exprs
  colnames(exprs) <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
  exprs <- asinh(exprs/5)
  exprs_set <- rbind(exprs_set, exprs)
  sample <- append(sample, rep(i, nrow(exprs)))
}
exprs_set$sample <- sample

# Set channels
pregating_channels <- c("DNA-1", "DNA-2")
lineage_channels <- c("CD45", "CD4", "CD20", "CD33", "CD123", "CD14", "IgM", "HLA-DR", "CD7", "CD3")
functional_channels <- c("pNFkB", "pp38", "pStat5", "pAkt", "pStat1", "pSHP2", "pZap70", "pStat3", "pSlp76", "pBtk", "pPlcg2", "pErk", "pLat", "pS6")
instrument_channels <- c("length", "Time", "sample", "BC1", "BC2", "BC3", "BC4", "BC5", "BC6", "BC7")

### Differential abundance of cells in clusters
# Sample 1000 random events form each sample
sub <- NULL
for(i in unique(exprs_set$sample)) {
  temp <- exprs_set[exprs_set$sample==i,]
  sub <- rbind(sub, temp[sample(nrow(temp), 1000, replace=FALSE),lineage_channels])
}
sample <- rep(unique(exprs_set$sample), each=1000)

# cluster
library(cytofkit)
clusters_pg <- cytof_cluster(xdata = sub, method = "Rphenograph")

# make frequency table
cell_counts <- matrix(nrow=max(clusters_pg), ncol=0)
for (i in unique(sample)) {
  sample_count <- table(clusters_pg[sample==i])
  column <- unname(sample_count[match(1:nrow(cell_counts), names(sample_count))])
  column[is.na(names(sample_count[match(1:nrow(cell_counts), names(sample_count))]))] <- 0
  cell_counts <- cbind(cell_counts, column)
}
conditions <- rep(1, 16)
conditions[grep("BCR-XL", fcs_files)] <- "BCR-XL"
conditions[grep("Reference", fcs_files)] <- "Reference"
colnames(cell_counts) <- conditions
# rownames(cell_counts) <- paste("cluster", 1:max(clusters_pg))
# cell_counts <- cell_counts[,order(colnames(cell_counts))]
# write.table(cell_counts, file="../DA_freq.txt", quote=FALSE, sep="\t")

# calculate differential abundance
library(edgeR)

dge <- DGEList(cell_counts, lib.size=rep(1000,16))
design <- model.matrix(~factor(conditions))
y <- estimateDisp(dge, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=2)

# visualize on t-SNE
library(Rtsne)
library(RColorBrewer)
library(scales)

set.seed(42)
tsne <- Rtsne(sub)

# par(mfrow=c(1,2))
layout(matrix(c(1,1,1,2,2,2,3), nrow=1,ncol=7))

getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
plot(tsne$Y, pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey", xlab="", ylab="", lwd=0.5, main="Phenograph clusters")

log2fc <- res$table$logFC
log2fc[res$table$PValue>0.05] <- 0
colpal <- colorRampPalette(c('red', "grey", 'blue'))
cols <- colpal(100)[as.numeric(cut(log2fc,breaks = 100))]
plot(tsne$Y, pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey", xlab="", ylab="", lwd=0.5, main="colored by log2FC")

legend_image <- as.raster(matrix(colpal(20), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = 'log2FC')
text(x=1.5, y = seq(0,1,l=3), labels = c(round(min(log2fc), digits = 2),0,round(max(log2fc),digits=2)))
rasterImage(legend_image, 0, 0, 1, 1)


# getPalette = colorRampPalette(brewer.pal(9, "Set1"))
# cols <- getPalette(max(clusters_pg))
# par(mfrow=c(1,2))
# plot(tsne$Y, pch=19, col="slategray4", xlab="", ylab="", lwd=0.5, main="Reference samples")
# points(tsne$Y[grepl("Reference", sample),1:2], pch=21, bg=alpha(col=cols[clusters_pg[grepl("Reference", sample)]],0.6), col="grey", xlab="", ylab="", lwd=0.5, cex=0.6)
# plot(tsne$Y, pch=19, col="slategray4", xlab="", ylab="", lwd=0.5, main="BCR-XL samples")
# points(tsne$Y[grepl("BCR-XL", sample),1:2], pch=21, bg=alpha(col=cols[clusters_pg[grepl("BCR-XL", sample)]],0.6), col="grey", xlab="", ylab="", lwd=0.5, cex=0.6)

### cydar
library(cydar)
library(ncdfFlow)

# Sample 1000 random events form each sample
sub <- list()
for(i in unique(exprs_set$sample)) {
  temp <- exprs_set[exprs_set$sample==i,]
  sub[[i]] <- temp[sample(nrow(temp), 1000, replace=FALSE),lineage_channels]
}
conditions <- rep(1, 16)
conditions[grep("BCR-XL", fcs_files)] <- "BCR-XL"
conditions[grep("Reference", fcs_files)] <- "Reference"
sample.id <- rep(c(1:16), each=1000)

# Counting cells into hyperspheres
cd_all <- prepareCellData(sub)
cd_all <- countCells(cd_all, tol=0.5)

# Testing for significant differences in abundance
library(edgeR)
y <- DGEList(assay(cd_all), lib.size=cd_all$totals)

keep <- aveLogCPM(y) >= aveLogCPM(5, mean(cd_all$totals))
cd <- cd_all[keep,]
y <- y[keep,]

design <- model.matrix(~factor(conditions))
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)
res <- glmQLFTest(fit, coef=2)

# Controlling the spatial FDR
qvals <- spatialFDR(intensities(cd), res$table$PValue)

is.sig <- qvals <= 0.05
summary(is.sig)

sig.coords <- intensities(cd)[is.sig,]
sig.res <- res$table[is.sig,]
coords <- prcomp(sig.coords)
plotCellLogFC(coords$x[,1], coords$x[,2], sig.res$logFC)

# Cluster and tsne + color by logFC
library(Rtsne)
library(cytofkit)
library(RColorBrewer)
library(scales)

set.seed(42)
tsne <- Rtsne(sub)

clusters_pg <- cytof_cluster(xdata = sub, method = "Rphenograph")

# color by cluster
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
plot(tsne$Y, pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey", xlab="", ylab="", lwd=0.5, main="Phenograph clusters")

# hypersphere assignment table
hslist <- data.frame(event=c(), hypersphere=c())
for(i in 1:length(cd_all@cellAssignments)) {
  hslist <- rbind(hslist, cbind(abs(cd_all@cellAssignments[[i]]), rep(i, length(cd_all@cellAssignments[[i]]))))
}
lapply(cd_all, function(x) cbind(unlist(x), rep(names(x)
hslist <- split(rep(length(cd@cellAssignments), lengths(cd@cellAssignments)), unlist(cd@cellAssignments))

# color by (significant) log2FC
res$table <- res$table[order(as.numeric(rownames(res$table))),]
log2fc <- res$table$logFC
log2fc[res$table$PValue>0.05] <- 0
colpal <- colorRampPalette(c('red','blue'))
cols <- colpal(100)[as.numeric(cut(LOG2FC,breaks = 100))]
# hypersphere assignments: cd@cellAssignments
# p val: res$table$PValue
# log2FC: res$table$logFC

# Check for duplicated events (remove if applicable - phenograph doesn't like duplicates)
table(duplicated(sub[,lineage_channels]))
sub <- sub[-which(duplicated(sub[,lineage_channels]))]
# Cluster using phenograph
clusters_pg <- cytof_cluster(xdata = sub[,lineage_channels], method = "Rphenograph")
# Make count table
cell_counts <- NULL
for(i in unique(sub$sample)) {
  table <- rep(0,max(clusters_pg)); names(table) <- 1:max(clusters_pg)
  table[names(table) %in% names(table(clusters_pg[sub$sample==i]))] <- unname(table(clusters_pg[sub$sample==i]))
  cell_counts <- cbind(cell_counts, table)
}
colnames(cell_counts) <- unique(sub$sample)
rownames(cell_counts) <- paste("cluster", 1:max(clusters_pg))

### Differential expression of proteins between clusters


## Cell type assignment

### MEM
# Load MEM library
library(MEM)
mem <- MEM(cbind(data[,colnames(data) %in% lm], cluster))
mem <- mem$MEM_matrix[[1]]

### flowCL
library("flowCL")

# Query (the output can be viewed using Res$Table)
Res <- flowCL("CD16-CD14+CD33+")


## Cellular hierarchies

### FlowSOM

# Calculate and visualize self organizing map with the parameters used by FlowSOM
library(kohonen)
som <- som(as.matrix(live[live$sample==fcs_files[1],lineage_channels]), grid=somgrid(xdim = 10, ydim = 10), dist.fcts="euclidean")

# Calculate codebook clusters
library(ConsensusClusterPlus)
nclust <- 15
cclust <- ConsensusClusterPlus(t(som$codes[[1]]), maxK = nclust, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(), plot = "pdf", verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = NULL)
codebook_clusters <- cclust[[nclust]]$consensusClass
codebook_list <- list()
for(i in 1:max(codebook_clusters)) {
  codebook_list[[i]] <- names(codebook_clusters[codebook_clusters==i])
}

# Create minimum spanning tree from codebook distance matrix
library(igraph)
graph  <- graph.adjacency(as.matrix(dist(som$codes[[1]])), weighted = TRUE)
mst <- mst(graph)

# Plot flowSOM-like MST output
layout <- layout.kamada.kawai(as.undirected(mst))
marker_palette = colorRampPalette(brewer.pal(9, "Set1"))
cluster_palette = colorRampPalette(brewer.pal(8, "Set2"))
marker_cols <- marker_palette(length(lineage_channels))
cluster_cols <- cluster_palette(max(codebook_clusters))
medians <- t(apply(som$codes[[1]], 1, function(x) x/sum(x)))
plot(as.undirected(mst), vertex.shape="star", vertex.label = NA, vertex.size = rep(10, 100), vertex.data = medians, vertex.cP = marker_palette(ncol(medians)), vertex.scale = TRUE, layout = layout, edge.lty = 1, mark.groups = codebook_list, mark.col = alpha(cluster_cols,0.4), mark.border = NA, mark.shape=1)
legend(x=1.2,y=1.1,legend=lineage_channels, col=marker_cols, pch=19, cex=0.7, bty="n", y.intersp=0.8)
# legend(x=1.2,y=1.1,legend=lineage_channels, col=marker_cols, pch=19, cex=1, bty="n", y.intersp=0.4)
legend(x=-1.8,y=1.1,legend=paste("cluster", 1:15), col=cluster_cols, pch=19, cex=0.7, bty="n", y.intersp=0.8)
# legend(x=-1.8,y=1.1,legend=paste("cluster", 1:15), col=cluster_cols, pch=19, cex=1, bty="n", y.intersp=0.4)

### SPADE
# Density-based downsampling: use dbscan
# A different approach is to use an unsupervised method, such as density-based clustering, This is relevant because the points we are really interested in are only the live, intact, singlet cells, which should be the dominating type of event in any sample, and tend to fall in regions of very high density. One such method is the so-called "Density-based spatial clustering of applications with noise" (DBSCAN), which is available for R in the package [dbscan](https://cran.r-project.org/web/packages/dbscan/). Both alternative approaches can be applied in a step-wise manner simulating the manual process, but potentially it should also work when using all the five dimensions at once (a problem with the speed may arise with this).

### FLOWSOM
###### SWITCH AT SOME POINT #######
## CHECK: https://www.r-bloggers.com/self-organising-maps-for-customer-segmentation-using-r/
## CLUSTERING IS JUST HIERACHICAL CLUSTERING ON CODEBOOK VECTORS USING ConsensusClusterPlus package and elbow method (unless cluster number is given)
library(kohonen)
som <- som(as.matrix(live[live$sample==fcs_files[1],lineage_channels]), grid=somgrid(xdim = 10, ydim = 10), dist.fcts="euclidean")
library(ConsensusClusterPlus)
nclust <- 15
cclust <- ConsensusClusterPlus(t(som$codes[[1]]), maxK = nclust, reps = 100, pItem = 0.9, pFeature = 1, title = tempdir(), plot = "pdf", verbose = FALSE, clusterAlg = "hc", distance = "euclidean", seed = NULL)
codebook_clusters <- cclust[[nclust]]$consensusClass
clusters_som <- unname(codebook_clusters[som$unit.classif])


#### CORRELATION MATRIX 
# Note that abundance is not the only useful metric. Cell functions overlap and generally interact a lot with each other, so a univariate metric is probably not *that* useful