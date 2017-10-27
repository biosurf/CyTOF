#!/usr/bin/env Rscript
##### "MATLAB" bead gating re-written in R #####
##### Christina Bligaard Pedersen, 2017 #####

# Installing packages
list.of.packages <- c("flowCore", "matrixStats", "pracma", "signal", "optparse", "e1071", "methods")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Loading libraries
loading = suppressWarnings(suppressMessages(lapply(list.of.packages, require, character.only = TRUE)))

# Sourcing functions
source("All_functions.R", chdir = TRUE)

# Commandline argument reader
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NA, 
              help="Directory containing input files", metavar="character"),
  make_option(c("-n", "--normalize"), type="character", default=NA, 
              help="Directory containing files to normalize input files to. If not specified, input folder will be used.", metavar="character"),
  #make_option(c("-v", "--version"), type="integer", default=3, 
  #            help="Set version of CyTOF used to generate data (1: v1, 2: v2, 3: Helios), default = 3.", metavar="integer"),
  make_option(c("-p", "--plot"), type="integer", default=1, 
              help="Generate plot showing bead gating? 1 = yes, 0 = no, default = 1.", metavar="integer")
)

# Parse command-line arguments
opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$input)) {
  stop("Command-line argument '-in' is required, use option -h to learn more.")
}


# Set version of CyTOF used to generate data (1: v1, 2: v2, 3: Helios) and co-factor. Currently version 3 is set automatically because of the automatic bead gating step.
v <- 3
if (v < 3) {c <- 5}; if (v == 3) {c <- 2}

# The script handles only DVS Beads (140,151,153,165,175), but may easily be extended
# Setting bead and DNA channels
channels = list(bead1 = 140, bead2 = 151, bead3 = 153, bead4 = 165, bead5 = 175, dna = 191)
n_beadchannels = 5

# Reading the input files
input <- get_file_list(opt$input)

# Set up output folders and log
log_file = generate_output(opt$input)
cat(c(format(Sys.time(), "%d-%h-%Y %H:%M:%S"),paste0("Normalizing files in folder ", opt$input)), sep = "\n",
              file=log_file)

# Identifying beads in all files
bead_indices = vector("list", input$num_files)
bead_data = matrix(list(), 1, input$num_files)
short_name = rep(NA, input$num_files)

# Loading bead SVM model
load("./Bead_SVM_model.Rdata")

# Identification of beads
for (i in 1:input$num_files) {
  file_data = read_file(input$filenames[i])
  short_name[i] = substr(basename(input$filenames[i]), 1, nchar(basename(input$filenames[i]))-4)
  file_bead_gating = SVM_bead_identify(file_data, channels, c, v, short_name[i], opt$input, opt$plot, bead_svm_model)
  bead_count = length(which(file_bead_gating$bead_gating == 'beads'))
  cat(c(basename(input$filenames[i]), file_bead_gating$output_string), sep = '\n', file=log_file, append=TRUE)

  if (bead_count == 0) {
    stop(paste0("No beads found in file ", input$filenames[i]), " - quitting the program.")
  }
  bead_indices[i] = list(which(file_bead_gating$bead_gating == 'beads'))
  bead_data[[1,i]] = file_data[file_bead_gating$bead_gating == 'beads',c("Time", colnames(file_data)[file_bead_gating$gating_channels[1:5]])]
}

all_beads = do.call(rbind, bead_data)


# Set baseline to which we're normalizing
# Setting folder to normalize against
if (is.na(opt$normalize)) {
  bead_means = colMedians(all_beads[,2:(n_beadchannels+1)])
} else {
  norm_files = get_file_list(opt$normalize)
  if (norm_files$num_files > 1){
    oldBead_vals = data.frame(matrix(ncol = 6, nrow = 0))
    
    for (i in 1:norm_files$num_files) {
      oldFile_data = read_file(norm_files$filenames[i])
      oldBead_vals = rbind(oldBead_vals, oldFile_data[,colnames(bead_data[[1]])])
    }
    
  } else {
    oldFile_data = read_file(norm_files$filenames)
    oldBead_vals = oldFile_data[,colnames(bead_data[[1]])]
  }
  bead_means = colMedians(as.matrix(oldBead_vals[,2:(n_beadchannels+1)]))
}


#### Time-dependent normalization ####
normedbeads = matrix(list(), 1, input$num_files)
smoothed = matrix(list(), 1, input$num_files)
normedsmoothed = matrix(list(), 1, input$num_files)

window_size = 200

for (i in 1:input$num_files) {
  file_data = read_file(input$filenames[i])
  file_bead_data = file_data[unlist(bead_indices[i]),c("Time", colnames(file_data)[file_bead_gating$gating_channels[1:5]])]
  
  # Smoothing beads
  smoothing = smooth_beads(window_size, file_bead_data)
  
  # Compute bead slopes
  bead_slopes = compute_bead_slopes(bead_means, smoothing$beads)
  
  # Checking that times are unique - in their code times are made unique (maybe this should be changed here)
  t = smoothing$time
  dt = diff(t)
  repeats = which(diff(smoothed$time)==0)
  if (length(repeats) > 0) {
    stop("Error - there are repeated times in the data")
  }
  
  # Correct channels (changes values of file_data)
  file_data = correct_channels(smoothing$time, bead_slopes, file_data)
  
  # Get smoothed beads and get normalized unsmoothed beads
  smoothed[[1,i]] = smoothing$beads
  normedbeads[[1,i]] = file_data[unlist(bead_indices[i]),c("Time", colnames(file_data)[file_bead_gating$gating_channels[1:5]])]
  
  # Smoothing normalized beads
  normedsmoothing = smooth_beads(window_size, normedbeads[[1,i]])
  normedsmoothed[[1,i]] = normedsmoothing$beads
  
  # Writing files
  ff = get_ff(input$filenames[i])
  ff2 = flowFrame(file_data[unlist(bead_indices[i]),], ff@parameters, ff@description)
  suppressWarnings(write.FCS(ff2, filename = paste0(opt$input, '/beads/', short_name[i], "_beads.fcs"), what="numeric", delimiter = "\\"))
  
  ff3 <- flowFrame(file_data, ff@parameters, ff@description)
  suppressWarnings(write.FCS(ff3, filename = paste0(opt$input, '/normed/', short_name[i], "_normalized.fcs"), what="numeric", delimiter = "\\"))
  
}


#### Plotting ####
# Plotting medians before and after
pdf(paste0(opt$input, '/normed/beads_before_and_after.pdf'), 16, 10)
par(cex.main = 2, cex.axis = 1.5, cex.lab = 1.5)
colors = c('blue', 'forestgreen', 'red', 'darkturquoise', 'magenta3')

s = matrix(0,nrow=input$num_files,ncol=1)
smoothedmeds = matrix(0, input$num_files, n_beadchannels)
normedmeds = matrix(0, input$num_files, n_beadchannels)

for (i in 1:input$num_files) {
  s[i,1] = nrow(smoothed[[i]])
  smoothedmeds[i,] = colMedians(bead_data[[i]][,2:(n_beadchannels+1)])
  normedmeds[i,] = colMedians(normedbeads[[i]][,2:(n_beadchannels+1)])
}

if (input$num_files > 1) {
  
  layout(matrix(c(1,1,1,1,1,2,3,3,3,3,3,4,3,3,3,3,3,5), nrow = 3, ncol = 6, byrow = TRUE))
  par(mar=c(5.1, 10.1, 4.1, 2.1))
  matplot(asinh(1/c*smoothedmeds), xlim = c(1-1/input$num_files,input$num_files+1/input$num_files), 
          ylim = c(min(asinh(1/c*smoothedmeds))-0.05, max(asinh(1/c*smoothedmeds))+0.05), type = 'o', 
          pch = 1, col = colors, lty=1, lwd=2, main = 'Before normalization', ylab = 'Transformed median intensities', xlab = '', xaxt='n')
  axis(1, labels = FALSE, tick = FALSE) 
  grid(nx = 0, ny = NULL)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "left",inset = 0,
         legend = colnames(bead_data[[i]][,2:(n_beadchannels+1)]), 
         col=colors, lwd=2, pch = 1, lty = 1, horiz = FALSE, y.intersp = 1.5, cex = 1.5)
  
  par(mar=c(30, 10.1, 4.1, 2.1))
  matplot(asinh(1/c*normedmeds), xlim = c(1-1/input$num_files,input$num_files+1/input$num_files), 
          ylim = c(min(asinh(1/c*smoothedmeds))-0.05, max(asinh(1/c*smoothedmeds))+0.05), type = 'o', 
          pch = 1, col = colors, lty=1, lwd=2, main = 'After normalization', ylab = 'Transformed median intensities', xlab = '', xaxt='n')
  par(mar=c(5.1, 10.1, 4.1, 2.1))
  axis(1, at = seq(1, input$num_files, 1), tick = FALSE, labels = FALSE) 
  text(x=seq(1, input$num_files, 1), y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]),
       labels=paste0(short_name, '.fcs'), srt=60, adj=1, xpd=TRUE, cex = 1.5)

  grid(nx = 0, ny = NULL)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "left",inset = 0,
         legend = colnames(bead_data[[i]][,2:(n_beadchannels+1)]), 
         col=colors, lwd=2, pch = 1, lty = 1, horiz = FALSE, y.intersp = 1.5, cex = 1.5)

} else {
  
  layout(matrix(c(1,1,1,1,1,2,3,3,3,3,3,4), nrow = 2, ncol = 6, byrow = TRUE))
  matplot(asinh(1/c*smoothed[[1]]), ylim = c(min(asinh(1/c*smoothed[[1]]))-0.05, max(asinh(1/c*smoothed[[1]]))+0.05), 
          type = 'l', pch = 19, col = colors, lty = 1, lwd = 2, main = paste0('Before normalization\n', 
          short_name[1], '.fcs'), ylab = 'Transformed intensities', xlab = 'Bead event number')
  
  grid(nx = 0, ny = NULL)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "left",inset = 0,
         legend = colnames(bead_data[[i]][,2:(n_beadchannels+1)]), 
         col=colors, lwd=2, horiz = FALSE, y.intersp = 1.5, cex = 1.5)
  
  matplot(asinh(1/c*normedsmoothed[[1]]), ylim = c(min(asinh(1/c*normedsmoothed[[1]]))-0.05, max(asinh(1/c*normedsmoothed[[1]]))+0.05), 
          type = 'l', pch = 19, col = colors, lty = 1, lwd = 2, main = paste0('After normalization\n', 
          short_name[1], '.fcs'), ylab = 'Transformed intensities', xlab = 'Bead event number')
  grid(nx = 0, ny = NULL)
  plot(1, type = "n", axes=FALSE, xlab="", ylab="")
  legend(x = "left",inset = 0,
         legend = colnames(bead_data[[i]][,2:(n_beadchannels+1)]), 
         col=colors, lwd=2, horiz = FALSE, y.intersp = 1.5, cex = 1.5)
}

invisible(dev.off())


# Calculating bead fractional ranges before/after normalization and printing to log
max_before = colMaxs(do.call(rbind, smoothed))
min_before = colMins(do.call(rbind, smoothed))

max_after = colMaxs(do.call(rbind, normedsmoothed))
min_after = colMins(do.call(rbind, normedsmoothed))

r_before=mean(max_before/min_before)
r_after=mean(max_after/min_after)
cat(c(paste0("Bead Fractional Range Before = ", round(r_before, 5)), paste0("Bead Fractional Range After = ", round(r_after, 5))), sep = '\n', file=log_file, append=TRUE)

close(log_file)

# Printing final message
cat(paste0("The script finished running for the files in the folder ", opt$input, ". Output files are found in subdirectories of ", opt$input, ".\n"))