######################################################
### The Flow Cytometry Standard (FCS) format {FCS} ###
######################################################

# Load the flowCore package
library(flowCore)

# Set working directory (directory where you placed the FCS files)
setwd("~/cytof_primer/fcs_files/")

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


################################
### CyTOF data preprocessing ###
################################

pregating_channels <- c("Bead", "DNA1", "DNA2", "Dead", "Event_length")
lineage_channels <- c("CD57", "CD19", "CD4", "CD8", "IgD", "CD11c", "CD16", "CD3", "CD38", "CD27", "CD14", "CXCR5", "CCR7", "CD45RA", "CD20", "CD127", "CD33", "CD28", "CD161", "TCRgd", "CD123", "CD56", "HLADR", "CD25")
instrument_channels <- c("Time", "Event_length", "Center", "Offset", "Width", "Residual", "sample")


########################################
### Dealing with "randomized" values ###
########################################

# Check if counts are randomized (expression matrix should only contain integers)
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
table(is.wholenumber(exprs[,c(pregating_channels, lineage_channels)]))

# You may revert back to original counts (only applicable to measured channels - i.e. exclude time, event length, etc.)
exprs <- cbind(ceiling(exprs[,c(pregating_channels, lineage_channels)]), exprs[,instrument_channels])


##############################
### ArcSinh transformation ###
##############################

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


##################
### Pre-gating ###
##################

# For the next step, set working directory to the normalized samples
setwd("~/cytof_primer/fcs_files/normed/")

# Read in your normalized files and concatenate in a data frame (in this example we subset each matrix to 10,000 random events)
fcs_files <- list.files(pattern='.fcs$', full=TRUE, ignore.case = TRUE)
exprs_set <- data.frame()
sample <- c()
for(i in fcs_files) {
  fcs <- read.FCS(filename=i, transformation=FALSE) 
  exprs <- fcs@exprs
  markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
  colnames(exprs)[which(!is.na(markers))] <- markers[which(!is.na(markers))]
  # exprs <- cbind(floor(exprs[,!colnames(exprs) %in% instrument_channels]), exprs[,colnames(exprs) %in% instrument_channels])
  exprs_set <- rbind(exprs_set, exprs[sample(nrow(exprs), 10000),])
  sample <- append(sample, rep(i, 10000))
}
exprs_set$sample <- sample
instrument_channels <- append(instrument_channels, "sample")

# Make arcsinh transformed expression matrix (excluding channels that should remain linear - e.g. time, event length, etc.)
cofac <- 2 # For v1 and v2 use a co-factor of 5, for Helios use a co-factor of 2
exprs_set_trans <-cbind(asinh(exprs_set[,c(pregating_channels, lineage_channels)]/cofac), exprs_set[,instrument_channels])


########################
### Batch correction ###
########################

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


########################################
### Gating for cells (beads vs. DNA) ###
########################################

# Load necessary packages
library(MASS)
library(RColorBrewer)

# Make color palette using color brewer
k <- 11; my.cols <- rev(brewer.pal(k, "RdYlBu"))

# Gate for cells (beads vs. DNA)
sub <- exprs_set_trans[exprs_set_trans$sample==fcs_files[4],c("Bead", "DNA1")]
z <- kde2d(x=sub[,1], y=sub[,2], n=50, h=max(bandwidth.nrd(sub[,1]), bandwidth.nrd(sub[,2])))

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 0; right <- 3; lower <- 6; upper <- 9
rect(left,lower,right,upper,border="red")
cells <- exprs_set_trans[exprs_set_trans[,"Bead"] < right,]; cells <- cells[cells[,"DNA1"] > lower & cells[,"DNA1"] < upper,]
paste("cells: ", round((nrow(cells)/nrow(exprs_set_trans))*100, digits=2), "%", sep="")


##############################################
### Gating for intact cells (DNA1 vs DNA2) ###
##############################################

# Gate for intact cells (DNA1 vs DNA2)
sub <- exprs_set_trans[,c("DNA1", "DNA2")]
z <- kde2d(sub[,1], sub[,2], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 6.3; right <- 8.2; lower <- 6.8; upper <- 8.85
rect(left,lower,right,upper,border="red")
intact <- cells[cells[,"DNA1"] > left & cells[,"DNA1"] < right,]; intact <- intact[intact[,"DNA2"] > lower & intact[,"DNA2"] < upper,]
paste("intact cells: ", round((nrow(intact)/nrow(cells))*100, digits=2), "%", sep="")


#################################################
### Gating for singlets (event length vs DNA) ###
#################################################

# Gate for singlets (event length vs. DNA)
sub <- intact[,c("Event_length", "DNA1")]
z <- kde2d(sub[,1], sub[,2], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 11; bottom <- min(sub[,2]); right <- 21; upper <- max(sub[,2])
rect(left,bottom,right,upper,border="red", lwd=3) # Draw cells gate (format: ll x,y; ur x,y)
singlets <- intact[intact[,"Event_length"] > left & intact[,"Event_length"] < right,]
paste("Intact singlets: ", round((nrow(singlets)/nrow(intact))*100, digits=2), "%", sep="")


######################################################
### Gating for live cells (live/dead stain vs DNA) ###
######################################################

# Gate for live
sub <- singlets[,c("Dead", "DNA1")]
z <- kde2d(sub[,1], sub[,2], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=colnames(sub)[1], ylab=colnames(sub)[2], xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 0; bottom <- min(sub[,2]); right <- 4; upper <- max(sub[,2])
rect(left,bottom,right,upper,border="red", lwd=3)
live <- singlets[singlets[,"Dead"] < right,]
paste("Live intact singlets: ", round((nrow(live)/nrow(singlets))*100, digits=2), "%", sep="")


for (i in fcs_files) {
  print(paste("Total recovery for ", i, ": ", round((nrow(live[live$sample==i,])/nrow(exprs_set_trans[exprs_set_trans$sample==i,]))*100, digits=2), "%", sep=""))
}

# Save FCS containing live intact singlets
ff2 <- flowFrame(xt_live_intact_singlets, ff@parameters, ff@description)
suppressWarnings(write.FCS(ff2, filename = paste("./",  file, "live_intact_singlets", sep=""), what="numeric", delimiter = "\\"))


####################################
### Principal component analysis ###
####################################

pca <- prcomp(live[live$sample==fcs_files[1],lineage_channels], scale. = TRUE)
plot(pca$x[,1:2], pch=21, bg=alpha("slategray4",0.6), col="grey")


####################
### t-SNE {tsne} ###
####################

table(duplicated(live[live$sample==fcs_files[1],lineage_channels]))
# No duplicates here, but if so, they can be removed:
live <- live[-which(duplicated(live[live$sample==fcs_files[1],lineage_channels]))]


library(Rtsne)
set.seed(42)
tsne <- Rtsne(live[live$sample==fcs_files[1],lineage_channels])
plot(tsne$Y[,1:2], pch=21, bg=alpha("slategray4",0.6), col="grey", xlab="", ylab="", lwd=0.5)


par(mfrow=c(2,2))
for(i in c("CD4", "CD8", "CD19", "CD56")) {
  cols_ramp <- colorRamp(c("slategray4", "red"))
  marker_exprs_scaled <- (live[live$sample==fcs_files[1],i] - min(live[live$sample==fcs_files[1],i]))/diff(range(live[live$sample==fcs_files[1],i]))
  cols <- rgb(cols_ramp(marker_exprs_scaled), maxColorValue = 256)
  plot(tsne$Y[,1:2], pch=21, bg=alpha(cols,0.6), col="grey", xlab="", ylab="", lwd=0.5, main=i)
}


##################
### Phenograph ###
##################

library(cytofkit)
clusters_pg <- cytof_cluster(xdata = live[live$sample==fcs_files[1],lineage_channels], method = "Rphenograph")


# Visualize results on PCA
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
plot(pca$x[,1:2], pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey")


# Visualize results on t-SNE plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_pg))
plot(tsne$Y[,1:2], pch=21, bg=alpha(col=cols[clusters_pg],0.6), col="grey", xlab="", ylab="", lwd=0.5)


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
  plot(tsne$Y[,1:2], pch=19, col="slategray4", xlab="", ylab="", lwd=0.5, main=paste("sample", i))
  points(tsne$Y[sub$sample==unique(sub$sample)[i],1:2], pch=21, bg=alpha(col=cols[clusters_pg[sub$sample==unique(sub$sample)[i]]],0.6), col="grey", xlab="", ylab="", lwd=0.5)
}


##########################
### FlowSOM clustering ###
##########################

library(cytofkit)
clusters_fc <- cytof_cluster(xdata = live[live$sample==fcs_files[1],lineage_channels], method = "FlowSOM", FlowSOM_k=25)


# Visualize results on t-SNE plot
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
cols <- getPalette(max(clusters_fc))
plot(tsne$Y[,1:2], pch=21, bg=alpha(col=cols[clusters_fc],0.6), col="grey", xlab="", ylab="", lwd=0.5)


################
### SCAFFoLD ###
################

# Loading required packages
require(devtools) || install.packages("devtools")
source("http://bioconductor.org/biocLite.R")
biocLite("flowCore")
library(devtools)
install_github("nolanlab/scaffold")
library(scaffold)

# Running SCAFFoLD (launching a GUI), please note that your current working directory should contain your files of interest before launching the GUI
scaffold.run()


######################
### Simple heatmap ###
######################

library(miscTools)
sub_matrix <- live[live$sample==fcs_files[1],lineage_channels]
cluster_matrix <- NULL
for(i in 1:max(clusters)) {
  cluster_matrix <- rbind(cluster_matrix, colMedians(sub_matrix[clusters==i,]))
}
library(gplots)
cols = colorRampPalette(c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[1]))
par(mar = c(2,2,2,2))
heatmap.2(t(cluster_matrix), col=cols, trace="none", density.info = "none", sepcolor = "white", sepwidth = c(0.001, 0.001), colsep=c(1:ncol(t(cluster_matrix))), rowsep=c(1:nrow(t(cluster_matrix))), xlab="cluster", ylab="channel")


# Load MEM library
library(MEM)
# Calculation of relative enrichment scores for each marker on each population using asinh transformation with co-factor 5
MEM_scores = MEM(clustered_data, transform=TRUE, cofactor=5, choose.ref=FALSE, IQR_thresh=NULL)
# Generation of heatmaps showing the marker profiles for each population - these are saved to files. 
build.heatmaps(MEM_scores, cluster.MEM="both", cluster.medians="none", display.thresh=0, output.files=TRUE)


##############
### Radviz ###
##############

library(Radviz)
cell.S <- make.S(colnames(sub_matrix))
cell.sim <- cosine(as.matrix(sub_matrix))
optim.cell <- do.optim(cell.S,cell.sim,iter=100,n=1000)
cell.S <- make.S(tail(optim.cell$best,1)[[1]])
cell.rv <- do.radviz(as.matrix(sub_matrix),cell.S)
layout(matrix(c(1,1,1,2), nrow=1,ncol=4))
plot(cell.rv, point.shape=19, point.color=alpha(cols[clusters], 0.6))
plot.new()
legend("center", legend=paste("cluster", 1:max(clusters)), col=cols, cex=1, pch=16, bty='n')


pop.rv <- do.radviz(cluster_matrix,cell.S)
size <- table(clusters)
layout(matrix(c(1,1,1,2), nrow=1,ncol=4))
bubbleRadviz(pop.rv, bubble.color=alpha(cols[1:nrow(cluster_matrix)],0.8), bubble.size=log(size[1:nrow(cluster_matrix)]), scale=0.2, decreasing=TRUE)
plot.new()
legend("center", legend=paste("cluster", 1:max(clusters)), col=cols, cex=1, pch=16, bty='n')


##############
### Citrus ###
##############

# Code modified from the original version by Robert Bruggner / bruggner@stanford.edu (output generated when using GUI)
# Compatible with Citrus version 0.8

rm(list = ls())
library("citrus",lib.loc=NULL)

# Use this line to limit the number of threads used by clustering
options("mc.cores"=1);

family = "classification"

# Change if you want to run from command line
# dataDirectory = "../"
dataDirectory = "/samples" 	 # The folder containing your files of interest
outputDirectory = file.path(dataDirectory,"citrusOutput")

clusteringColumns = c("CD3(110:114)Dd","CD45(In115)Dd","CD4(Nd145)Dd","CD20(Sm147)Dd","CD33(Nd148)Dd","CD123(Eu151)Dd","CD14(Gd160)Dd","IgM(Yb171)Dd","HLA-DR(Yb174)Dd","CD7(Yb176)Dd")		# Parameters to use in the clustering step - the channels used for manual gating
medianColumns=c("pNFkB(Nd142)Dd","pp38(Nd144)Dd","pStat5(Nd150)Dd","pAkt(Sm152)Dd","pStat1(Eu153)Dd","pSHP2(Sm154)Dd","pZap70(Gd156)Dd","pStat3(Gd158)Dd","pSlp76(Dy164)Dd","pBtk(Er166)Dd","pPlcg2(Er167)Dd","pErk(Er168)Dd","pLat(Er170)Dd","pS6(Yb172)Dd")	
transformColumns = c("CD3(110:114)Dd","CD45(In115)Dd","BC1(La139)Dd","BC2(Pr141)Dd","pNFkB(Nd142)Dd","pp38(Nd144)Dd","CD4(Nd145)Dd","BC3(Nd146)Dd","CD20(Sm147)Dd","CD33(Nd148)Dd","pStat5(Nd150)Dd","CD123(Eu151)Dd","pAkt(Sm152)Dd","pStat1(Eu153)Dd","pSHP2(Sm154)Dd","pZap70(Gd156)Dd","pStat3(Gd158)Dd","BC4(Tb159)Dd","CD14(Gd160)Dd","pSlp76(Dy164)Dd","BC5(Ho165)Dd","pBtk(Er166)Dd","pPlcg2(Er167)Dd","pErk(Er168)Dd","BC6(Tm169)Dd","pLat(Er170)Dd","IgM(Yb171)Dd","pS6(Yb172)Dd","HLA-DR(Yb174)Dd","BC7(Lu175)Dd","CD7(Yb176)Dd","DNA-1(Ir191)Dd","DNA-2(Ir193)Dd")			#  Markers to transform
transformCofactor = 5 		# Normal choice for mass cytometry data
scaleColumns = c("CD3(110:114)Dd","CD45(In115)Dd","BC1(La139)Dd","BC2(Pr141)Dd","pNFkB(Nd142)Dd","pp38(Nd144)Dd","CD4(Nd145)Dd","BC3(Nd146)Dd","CD20(Sm147)Dd","CD33(Nd148)Dd","pStat5(Nd150)Dd","CD123(Eu151)Dd","pAkt(Sm152)Dd","pStat1(Eu153)Dd","pSHP2(Sm154)Dd","pZap70(Gd156)Dd","pStat3(Gd158)Dd","BC4(Tb159)Dd","CD14(Gd160)Dd","pSlp76(Dy164)Dd","BC5(Ho165)Dd","pBtk(Er166)Dd","pPlcg2(Er167)Dd","pErk(Er168)Dd","BC6(Tm169)Dd","pLat(Er170)Dd","IgM(Yb171)Dd","pS6(Yb172)Dd","HLA-DR(Yb174)Dd","BC7(Lu175)Dd","CD7(Yb176)Dd","DNA-1(Ir191)Dd","DNA-2(Ir193)Dd")			# @Lars, help I don't know what this is...

minimumClusterSizePercent = 0.05 	# Ignore clusters with less than fileSampleSize* minimumClusterSizePercent events
modelTypes = c("pamr")	 		   # Association model selection (see description below, options: "pamr", "glmnet", "sam")
fileSampleSize = 10000		       # How many events to use from each sample
nFolds = 4                           # Cross-validation folds to use
featureType = c("medians") 	 	  # How to characterize the clusters (options: "medians", "abundances")

# Setting files to run on and grouping them (Group 1 is the reference, Group 2 is the stimulated data)
fileList = data.frame(defaultCondition=c("PBMC8_30min_patient1_Reference.fcs","PBMC8_30min_patient2_Reference.fcs","PBMC8_30min_patient3_Reference.fcs","PBMC8_30min_patient4_Reference.fcs","PBMC8_30min_patient5_Reference.fcs","PBMC8_30min_patient6_Reference.fcs","PBMC8_30min_patient7_Reference.fcs","PBMC8_30min_patient8_Reference.fcs","PBMC8_30min_patient1_BCR-XL.fcs","PBMC8_30min_patient2_BCR-XL.fcs","PBMC8_30min_patient3_BCR-XL.fcs","PBMC8_30min_patient4_BCR-XL.fcs","PBMC8_30min_patient5_BCR-XL.fcs","PBMC8_30min_patient6_BCR-XL.fcs","PBMC8_30min_patient7_BCR-XL.fcs","PBMC8_30min_patient8_BCR-XL.fcs"))

labels = as.factor(c("Group 1","Group 1","Group 1","Group 1","Group 1","Group 1","Group 1","Group 1","Group 2","Group 2","Group 2","Group 2","Group 2","Group 2","Group 2","Group 2"))

# Running citrus
results = citrus.full(
            fileList=fileList,
            labels=labels,
            clusteringColumns=clusteringColumns,
            dataDirectory=dataDirectory,
            outputDirectory=outputDirectory,
            family=family,
            modelTypes=modelTypes,
            nFolds=nFolds,
            
            fileSampleSize=fileSampleSize,
            featureType=featureType,
            minimumClusterSizePercent=minimumClusterSizePercent,
            transformColumns=transformColumns,
            transformCofactor=transformCofactor,
            scaleColumns=scaleColumns,
            medianColumns=medianColumns
)

# Plotting results
plot(results,outputDirectory)


##############
### flowCL ###
##############

# Loading packages
source("http://bioconductor.org/biocLite.R")
biocLite("flowCL")
library("flowCL")

# Query (the output can be viewed using Res$Table)
Res <- flowCL("CD16-CD14+CD33+")


###############
### FlowSOM ###
###############

# Calculate and visualize self organizing map with the parameters used by FlowSOM
library(kohonen)
som <- som(as.matrix(live[live$sample==fcs_files[1],lineage]), grid=somgrid(xdim = 10, ydim = 10), dist.fcts="euclidean")
plot(som, type="codes", codeRendering="segments")

# Minimum spanning tree visualization
library(igraph)
graph  <- graph.adjacency(as.matrix(dist(som$codes[[1]])), weighted = TRUE)
mst <- mst(graph)
plot(as.undirected(mst), vertex.size=10, mode="undirected")
