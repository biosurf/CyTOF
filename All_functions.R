#### Getting list of input files ####
get_file_list = function(data_dir) {
  files <- list.files(path = data_dir, pattern='.fcs', ignore.case=TRUE, full=TRUE)
  
  # Error check
  if (length(files) == 0) {
    stop("You have specified a folder with no fcs-files.")
  } else if (nrow(file.info(files)[file.info(files)$size == 0,]) != 0) {
    stop("Empty fcs-files are present in the directory - remove these and try again.")
  } else {
    num_files = length(files)
  }
  
  return(list(filenames = files, num_files = num_files))
}


#### Setting up output directories and log ####
generate_output = function(data_dir) {
  dir.create(file.path(data_dir, "beads"), showWarnings = FALSE)
  dir.create(file.path(data_dir, "normed"), showWarnings = FALSE)
  log_file = file(paste0(data_dir, '/normed/normalizer.log'), 'w')
  
  return(log_file)
}


#### Reading fcs file and getting the data matrix ####
read_file = function(filename) {
  ff = read.FCS(filename, transformation=FALSE)
  markers = as.vector(ff@parameters@data[,2])
  names(markers) = as.vector(ff@parameters@data[,1])
  file_data = ff@exprs
  
  return(file_data)
}


#### Identification of beads ####
SVM_bead_identify = function(file_data, channels, c, v, short_name, directory, plot, bead_svm_model) {

  # Transforming data
  trans_data = file_data[,colnames(file_data) %in% c("Time", "Event_length")]
  trans_data = cbind(trans_data, asinh(file_data[,!colnames(file_data) %in% c("Time", "Event_length")]/c))
  if (v < 3) {file_data[file_data<0] <- 0}
  
  # Predicting bead events (using 'normalized' data)
  trans_data_norm <- cbind(trans_data[,1:2], apply(trans_data[,3:ncol(trans_data)], 2, function(x)(x-min(x))/(max(x)-min(x))))
  gating_channels <- c(grep(channels$bead1, colnames(trans_data)), grep(channels$bead2, colnames(trans_data)), grep(channels$bead3, colnames(trans_data)), grep(channels$bead4, colnames(trans_data)), grep(channels$bead5, colnames(trans_data)), grep(channels$dna, colnames(trans_data)))
  bead_gating <- predict(bead_svm_model, trans_data_norm[,gating_channels])
  # Plotting bead gating
  output_string = plot_beads(trans_data, gating_channels, bead_gating, short_name, directory, plot)
  
  return(list(bead_gating = bead_gating, gating_channels = gating_channels, output_string = output_string))
}


#### Plotting bead gating ####
plot_beads = function(data, channels, gating, filename, directory, plot) {
  output_string = ""
  
  if (plot == 1) {
    x <- channels[1]
    y <- channels[2]
    z <- channels[3]
    v <- channels[4]
    w <- channels[5]
    q <- channels[6]
    
    pdf(paste0(directory, '/beads/', filename, '.pdf'), width = 20, height = 5)
    par(mfrow = c(1, 5))
    plot(data[,c(x, q)], pch=".", cex=0.1, col="gray", xlab=colnames(data)[x], ylab=colnames(data)[q], xlim=c(0,max(data[,c(x, q)][,1])), ylim=c(0,max(data[,c(x, q)][,2])), main='Bead gate 1')
    points(data[gating == 'beads',c(x, q)], pch=".", cex=0.1, col="red")
    
    plot(data[,c(y, q)], pch=".", cex=0.1, col="gray", xlab=colnames(data)[y], ylab=colnames(data)[q], xlim=c(0,max(data[,c(y, q)][,1])), ylim=c(0,max(data[,c(y, q)][,2])), main='Bead gate 2')
    points(data[gating == 'beads',c(y, q)], pch=".", cex=0.1, col="red")
    
    plot(data[,c(z, q)], pch=".", cex=0.1, col="gray", xlab=colnames(data)[z], ylab=colnames(data)[q], xlim=c(0,max(data[,c(z, q)][,1])), ylim=c(0,max(data[,c(z, q)][,2])), main='Bead gate 3')
    points(data[gating == 'beads',c(z, q)], pch=".", cex=0.1, col="red")
    
    plot(data[,c(v, q)], pch=".", cex=0.1, col="gray", xlab=colnames(data)[v], ylab=colnames(data)[q], xlim=c(0,max(data[,c(v, q)][,1])), ylim=c(0,max(data[,c(v, q)][,2])), main='Bead gate 4')
    points(data[gating == 'beads',c(v, q)], pch=".", cex=0.1, col="red")
    
    plot(data[,c(w, q)], pch=".", cex=0.1, col="gray", xlab=colnames(data)[w], ylab=colnames(data)[q], xlim=c(0,max(data[,c(w, q)][,1])), ylim=c(0,max(data[,c(w, q)][,2])), main='Bead gate 5')
    points(data[gating == 'beads',c(w, q)], pch=".", cex=0.1, col="red")
    dev.off()
    
    output_string = paste0("   A plot of the gated beads was saved as ", directory, '/beads/', filename, ".pdf.\n")
  }
  
  n_beads = nrow(data[gating == 'beads',])
  p_beads = n_beads/(nrow(data))*100
  output_string = paste0(output_string, "   ", n_beads, " beads found (", round(p_beads, 2), "% of all events)", ".")
  
  return(output_string)
}


#### Smoothing beads in windows ####
smooth_beads = function(window_size, file_bead_data) {
  
  num_events = nrow(file_bead_data)
  
  if (num_events >= window_size) {
    smoothed_beads = matrix(0,nrow=(num_events-window_size),ncol=n_beadchannels)
    smoothed_time = matrix(0,nrow=(num_events-window_size),ncol=1)
    
    for (i in 1:(num_events-window_size)) {
      smoothed_beads[i,] = colMedians(file_bead_data[i:(i+window_size-1),2:(n_beadchannels+1)])   # Here, the original code does not subtract 1 from the window size. I do it because it is supposed to be a 200-bead window (not a 201-bead window). It probably does not make a large difference.
      smoothed_time[i] = median(file_bead_data[i:(i+window_size-1),1])
    }

  } else {
    smoothed_beads = colMedians(file_bead_data[,2:(n_beadchannels+1)])
    smoothed_time = median(file_bead_data[,1])
  }
  
  return(list(beads = smoothed_beads, time = smoothed_time))
}


#### Computation of bead slopes in a window ####
compute_bead_slopes = function(baseline, smoothed_beads) {
  s = nrow(smoothed_beads)
  y = repmat(baseline, s, 1)
  
  bead_slopes = rowSums(smoothed_beads*y) / rowSums(smoothed_beads^2)
  
  return(bead_slopes)
}


#### Correct channels according to bead slopes ####
correct_channels = function(smoothed_time, bead_slopes, file_data) {
  if (length(smoothed_time) > 1) {
    ir = interp1(smoothed_time[,1], bead_slopes, file_data[,"Time"], 'linear')
    
    # Find border time points and extrapolate
    tstart2 = head(which(!is.na(ir)), 1)
    ir[1:tstart2] = ir[tstart2]
    
    tend = tail(which(!is.na(ir)), 1)
    ir[tend:length(ir)] = ir[tend]
  } else {
    ir = bead_slopes * matrix(1,nrow=nrow(file_data),ncol=1)
  }
  
  for (j in which(!colnames(file_data) %in% c("Time", "Event_length"))) {
    file_data[,j] = file_data[,j]*ir   # Performing the normalization
  }

  return(file_data)
}


#### Reading an FCS file and getting information ####
get_ff = function(filename) {
  ff = read.FCS(filename, transformation=FALSE)
  return(ff)
}
