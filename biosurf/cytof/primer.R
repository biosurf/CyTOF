## The Flow Cytometry Standard (FCS) format

# Load the flowCore package
library(flowCore)

# Set working directory (directory where you placed the FCS files)
setwd("~/Dropbox/Research/Papers/Ongoing/CyTOF primer/CyTOF_wiki/data/raw/")

# Read an fcs file
fcs_files <- list.files(pattern='.fcs$', full=TRUE, ignore.case = TRUE)
fcs <- read.FCS(filename=fcs_files[1], transformation=FALSE)

# Extract expression matrix
exprs <- fcs@exprs

# Make channels human readable using information in the parameter data slot (note: formats may vary - examine fcs@parameters@data)
markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs@parameters@data$desc))
colnames(exprs)[which(!is.na(markers))] <- temp[which(!is.na(markers))]

# Merging fcs files
source("~/Dropbox/Git_repos/r-cytof/concatenateDirectoryFiles.R")
dir.create("./combined")
cytofCore.concatenateDirectoryFiles(inputDir="./",outputDir="./combined",pattern=NULL,overwrite=F,timeParam="time")

## Data transformations

### Randomization

# Check if counts are randomized (whole numbers are )
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
table(is.wholenumber(exprs[,!colnames(exprs) %in% c("Time", "Event_length", "Center", "Offset", "Width", "Residual")]))

# Revert back to original count (only applicable to measured channels)
exprs <- cbind(exprs, floor(exprs[,!colnames(exprs) %in% c("Time", "Event_length", "Center", "Offset", "Width", "Residual")]), exprs[,colnames(exprs) %in% c("Time", "Event_length", "Center", "Offset", "Width", "Residual")])

### ArcSinh transformation

# Set version of CyTOF used to generate data (1: v1, 2:v2, 3: Helios) and the co-factor

v <- 3
if (v < 3) {c <- 5}; if (v == 3) {c <- 2}

# Make arcsinh transformed expression matrix (excluding channels that should remain linear - e.g. Time and Event_length)

exprs_trans <- exprs[,colnames(exprs) %in% c("Time", "Event_length", "Center", "Offset", "Width", "Residual")]
exprs_trans <- cbind(exprs_trans, asinh(exprs[,!colnames(exprs) %in% c("Time", "Event_length", "Center", "Offset", "Width", "Residual")]/c))
if (v < 3) {xt[xt<0] <- 0}
