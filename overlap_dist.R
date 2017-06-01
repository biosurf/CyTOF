# n=3
# s=spectra
# l=lasers
# plot=TRUE

overlap <- function(n, s, l, r, plot) {
  if (r == "b") {
    # for (n in 2:16) {
    suppressWarnings(library(matrixStats))
    source("~/Dropbox/Git_repos/biosurf/flow/wavelength_to_rgb.R")
    
    auc_spectra <- list()
    for(i in 1:(length(s))) {
      auc_spectra[[i]] <- sum(s[[i]]$Emission)
    }
    
    peak_dist <- matrix(ncol=length(s), nrow=0)
    for (i in 1:(length(s))) {
      row <- c()
      for (ii in 1:(length(s))) {
        dist <- sqrt(((s[[ii]]$Wavelength[which.max(s[[ii]]$Emission)] - s[[i]]$Wavelength[which.max(s[[i]]$Emission)])^2) + (((max(s[[ii]]$Emission) - max(s[[i]]$Emission))^2)^2))
        row <- append(row, dist)
      }
      peak_dist <- rbind(peak_dist, as.numeric(row))
    }
    peak_dist[lower.tri(peak_dist, diag = TRUE)] <- 0
    
    # Pick first two fluorophores and find min and max
    max_dist_pair <- as.vector(which(peak_dist == max(peak_dist), arr.ind = TRUE))
    min <- s[[max_dist_pair[1]]]$Wavelength[which.max(s[[max_dist_pair[1]]]$Emission)]
    max <- s[[max_dist_pair[2]]]$Wavelength[which.max(s[[max_dist_pair[2]]]$Emission)]
    
    # For each n above 2, find eqidistant spectra peaks
    peak_wl <- c()
    for(i in 1:length(s)) {
      peak_wl <- append(peak_wl, s[[i]]$Wavelength[which.max(s[[i]]$Emission)])
    }
    best_combi <- max_dist_pair
    if(n>2) {
      max_dist <- max-min
      eqi_points <- seq(min+(max-min)/(n-1),max-(max-min)/(n-1),(max-min)/(n-1))
      for(i in 1:(n-2)) {
        # Find peak closest to each eqidistant points
        best_combi <- append(best_combi, which.min(abs(eqi_points[i]-peak_wl)))
      }
    }
    
    # Plot results
    # pdf(file = paste("~/Dropbox/Research/Projects/Ongoing/spectral_overlapper/n=", n, "_simple.pdf", sep=""), height=8, width=8)
    if(plot==TRUE) {
      plot(1, type="n", xlab="Wavelength (nm)", ylab="Relative intensity (%)", xlim=c(300, 900), ylim=c(0, 100), yaxs="i", xaxs="i")
      grid(14,9, lty = 1, lwd = 0.5)
      for(i in 1:length(best_combi)) {
        points(s[[best_combi[i]]]$Wavelength, s[[best_combi[i]]]$Excitation, col=wavelength_to_rgb(s[[best_combi[i]]]$Wavelength[which.max(s[[best_combi[i]]]$Emission)]), type="l", lty=2)
      }
      for(i in 1:length(best_combi)) {
        polygon(s[[best_combi[i]]]$Wavelength, s[[best_combi[i]]]$Emission, col=adjustcolor(col=wavelength_to_rgb(s[[best_combi[i]]]$Wavelength[which.max(s[[best_combi[i]]]$Emission)]),alpha.f = 0.5), border = adjustcolor(col=wavelength_to_rgb(s[[best_combi[i]]]$Wavelength[which.max(s[[best_combi[i]]]$Emission)]),alpha.f = 2))
      }
      for(i in 1:length(best_combi)) {
        text(s[[best_combi[i]]]$Wavelength[which.max(s[[best_combi[i]]]$Emission)], 20, labels = names(s)[best_combi[i]], cex = 0.7, srt=90)
      }
      for(i in 1:length(l)) {
        abline(v=l[i], lty=2, col="dark grey")
        text(l[i], 95, labels=paste(l[i], "nm"), cex=0.5, srt=45)
      }
    }
    # dev.off()
    # }
  }
  if (r == "b") {
    ## WRITE ADVANCED OPTION
  }
}
