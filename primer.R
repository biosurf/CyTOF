################################################
### The Flow Cytometry Standard (FCS) format ###
################################################

# Load the flowCore package
library(flowCore)

# Set working directory (directory where you placed the FCS files)
setwd("~/Dropbox/Research/Papers/Ongoing/CyTOF primer/CyTOF_wiki/data/raw/")

# Read an fcs file
fcs_files <- list.files(pattern='.fcs$', full=TRUE, ignore.case = TRUE)
fcs <- read.FCS(filename=fcs_files[1], transformation=FALSE)

# Extract expression matrix
exprs <- fcs@exprs

# Make channels human readable using information in the parameter data slot (fcs@parameters@data)
markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
colnames(exprs)[which(!is.na(markers))] <- temp[which(!is.na(markers))]

# Merging fcs files
source("~/Dropbox/Git_repos/r-cytof/concatenateDirectoryFiles.R")
dir.create("./combined")
cytofCore.concatenateDirectoryFiles(inputDir="./",outputDir="./combined",pattern=NULL,overwrite=F,timeParam="time")


#####################
### Randomization ###
#####################

# Check if counts are randomized

# Revert back to original count



##############################
### ArcSinh transformation ###
##############################

# Set version of CyTOF used to generate data (1: v1, 2:v2, 3: Helios) and the co-factor
v <- 3
if (v < 3) {c <- 5}; if (v == 3) {c <- 2}

# Make arcsinh transformed expression matrix (with the exception of time and event_length, which should remain linear)
xt <- x[,colnames(x) %in% c("Time", "Event_length")]
xt <- cbind(xt, asinh(x[,!colnames(x) %in% c("Time", "Event_length")]/c))
if (v < 3) {xt[xt<0] <- 0}



# Transformation by channel
xx <- melt(expr_count[1:2000,])
colnames(xx) <- c("Marker", "Intensity")
ggplot(xx, aes(x=Marker, y=Intensity)) + geom_boxplot() + ggtitle("Raw")

xx <- melt(asinh(expr_count[1:2000,]/5))
colnames(xx) <- c("Marker", "Intensity")
ggplot(xx, aes(x=Marker, y=Intensity)) + geom_boxplot() + ggtitle("ArcSinh transformed")


##################
### Pre-gating ###
##################

# Set the channels you used for beads, DNA1, DNA2, and live/dead (cisplatin)
bead <- 140; dna1 <- 191; dna2 <- 193; dead <- 115

# Gate for cells
x <- names(markers)[grep(bead, names(markers), ignore.case = TRUE)]
y <- names(markers)[grep(dna1, names(markers), ignore.case = TRUE)]
sub <- xt[,c(x, y)]

z <- kde2d(xt[,x], xt[,y], n=50)  ### DET HER VIRKER IKKE!!! 
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)  ### DET HER VIRKER IKKE!!! 

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=paste(unname(markers[x]), " (", names(markers[x]),")", sep=""), ylab=paste(unname(markers[y]), " (", names(markers[y]),")", sep=""), xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,2])))

# Set gate boundaries - the values should be manually adjusted
left <- 0; right <- 3; lower <- 6; upper <- 9
rect(left,lower,right,upper,border="red")
xt_cells <- xt[xt[,x] < right,]; xt_cells <- xt_cells[xt_cells[,y] > lower & xt_cells[,y] < upper,]
paste("cells: ", round((nrow(xt_cells)/nrow(xt))*100, digits=2), "%", sep="")
text((left+(right-left)/2),(lower+(upper-lower)/2)-2, paste("cells: ", round((nrow(xt_cells)/nrow(xt))*100, digits=2), "%", sep=""), col="red")



# Gate for intact cells
x <- names(markers)[grep(dna2, names(markers), ignore.case = TRUE)]
y <- names(markers)[grep(dna1, names(markers), ignore.case = TRUE)]
sub <- xt_cells[,c(x, y)]
z <- kde2d(xt_cells[,x], xt_cells[,y], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=paste(unname(markers[x]), " (", names(markers[x]),")", sep=""), ylab=paste(unname(markers[y]), " (", names(markers[y]),")", sep=""), xlim=c(0,max(sub[,1])), ylim=c(0,max(sub[,1])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 7.25; right <- 8.4; lower <- 6.6; upper <- 7.85
rect(left,lower,right,upper,border="red")
xt_intact <- xt_cells[xt_cells[,x] > left & xt_cells[,x] < right,]; xt_intact <- xt_intact[xt_intact[,y] > lower & xt_intact[,y] < upper,]
paste("intact cells: ", round((nrow(xt_intact)/nrow(xt_cells))*100, digits=2), "%", sep="")
text((left+(right-left)/2),(lower+(upper-lower)/2)-2, paste("intact cells: ", round((length(xt_intact)/length(xt_cells))*100, digits=2), "%", sep=""), col="red")



# Gate for singlets
x <- names(markers)[grep("length", names(markers), ignore.case = TRUE)]
y <- names(markers)[grep(dna1, names(markers), ignore.case = TRUE)]
sub <- xt_intact[,c(x, y)]
z <- kde2d(xt_intact[,x], xt_intact[,y], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=x, ylab=paste(unname(markers[y]), " (", names(markers[y]),")", sep=""), xlim = c(0,max(xt[,x])), ylim = c(0,max(xt[,y])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 12; right <- 23
rect(left,lower,right,upper,border="red")
xt_intact_singlets <- xt_intact[xt_intact[,x] >= left & xt_intact[,x] <= right,]
paste("Intact singlets: ", round((nrow(xt_intact_singlets)/nrow(xt_intact))*100, digits=2), "%", sep="")
text((left+(right-left)/2),(lower+(upper-lower)/2-1),paste("intact singlets: ", round((length(xt_intact_singlets)/length(xt_intact))*100, digits=2), "%", sep=""), col="red")



# Gate for live
x <- names(markers)[grep(dead, names(markers), ignore.case = TRUE)]
y <- names(markers)[grep(dna1, names(markers), ignore.case = TRUE)]
sub <- xt_intact_singlets[,c(x, y)]
z <- kde2d(xt_intact_singlets[,x], xt_intact_singlets[,y], n=50)

# Generating a plot for gating visualization
plot(sub, pch=".", cex=0.1, col="gray", xlab=paste(unname(markers[x]), " (", names(markers[x]),")", sep=""), ylab=paste(unname(markers[y]), " (", names(markers[y]),")", sep=""), xlim = c(0,max(xt[,c(x)])), ylim = c(0,max(xt[,c(y)])))
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)

# Set gate boundaries - the values should be manually adjusted
left <- 0; right <- 4.6
rect(left,lower,right,upper,border="red") # Draw cells gate (format: ll x,y; ur x,y)
xt_live_intact_singlets <- xt_intact_singlets[xt_intact_singlets[,x] < right,]
paste("Live intact singlets: ", round((nrow(xt_live_intact_singlets)/nrow(xt_intact_singlets))*100, digits=2), "%", sep="")
text((left+(right-left)/2),(lower+(upper-lower)/2-1),paste("live intact singlets: ", round((length(xt_live_intact_singlets)/length(xt_intact_singlets))*100, digits=2), "%", sep=""), col="red")



paste("Total recovery: ", round((nrow(xt_live_intact_singlets)/nrow(xt))*100, digits=2), "%", sep="")

# Save FCS containing live intact singlets
ff2 <- flowFrame(xt_live_intact_singlets, ff@parameters, ff@description)
suppressWarnings(write.FCS(ff2, filename = paste("./",  file, "live_intact_singlets", sep=""), what="numeric", delimiter = "\\"))


##############
### Radviz ###
##############

expr_norm <- apply(expr_bc[!is.na(pd_new$cluster),3:12],2,do.L,fun=function(x) quantile(x,c(0.025,0.975)))
cell.S <- make.S(dimnames(expr_norm)[[2]])
cell.sim <- cosine(expr_norm)
optim.cell <- do.optim(cell.S,cell.sim,iter=100,n=1000)
cell.S <- make.S(tail(optim.cell$best,1)[[1]])
cell.rv <- do.radviz(expr_norm,cell.S)
plot(cell.rv, point.shape=19, point.color=cols[pd_new$cluster[!is.na(pd_new$cluster)]])


###########
### MEM ###
###########

# Load MEM library
library(MEM)
# Calculation of relative enrichment scores for each marker on each population using asinh transformation with co-factor 5
MEM_scores = MEM(clustered_data, transform=TRUE, cofactor=5, choose.ref=FALSE, IQR_thresh=NULL)
# Generation of heatmaps showing the marker profiles for each population - these are saved to files. 
build.heatmaps(MEM_scores, cluster.MEM="both", cluster.medians="none", display.thresh=0, output.files=TRUE)
