# n=6
# s=spectra
# l=lasers
# plot=TRUE

overlap <- function(n, s, l, r, plot) {
  if (r == "b") {
    suppressWarnings(library(matrixStats))
    source("~/Dropbox/Git_repos/flow/wavelength_to_rgb.R")
    
    # Calculate auc for each fluorophore [CAN BE PRE-CALCULATED]
    auc_spectra <- list()
    for(i in 1:(length(s))) {
      auc_spectra[[i]] <- sum(s[[i]]$Emission)
    }
    
    # Calculate auc_overlaps between fluorophores [CAN BE PRE-CALCULATED]
    auc_overlaps <- matrix(ncol=length(s), nrow=0)
    for (i in 1:(length(s))) {
      row <- c()
      for (ii in 1:(length(s))) {
        # Wave lengths where i and ii overlap
        wl_ol <- intersect(s[[i]]$Wavelength[which(s[[i]]$Emission>0)], s[[ii]]$Wavelength[which(s[[ii]]$Emission>0)]); if(length(wl_ol) == 0) {wl_ol <- 0} else {wl_ol <- c(min(wl_ol)-1, wl_ol, max(wl_ol)+1)}
        wl_ol_vec <- which(s[[i]]$Wavelength %in% wl_ol)
        # Loss of i
        loss_i <- s[[i]]$Emission[wl_ol_vec] - s[[ii]]$Emission[wl_ol_vec]
        loss_i <- sum(s[[i]]$Emission[wl_ol_vec[which(loss_i <= 0)]])
        loss_i <- loss_i/auc_spectra[[i]] 
        # Loss of ii
        loss_ii <- s[[ii]]$Emission[wl_ol_vec] - s[[i]]$Emission[wl_ol_vec]
        loss_ii <- sum(s[[ii]]$Emission[wl_ol_vec[which(loss_ii <= 0)]])
        loss_ii <- loss_ii/auc_spectra[[ii]] 
        # Append combined fractional loss to row
        row <- append(row, loss_i+loss_ii)
      }
      auc_overlaps <- rbind(auc_overlaps, as.numeric(row))
    }
    
    # List combinations 
      # Can be pre-calculated
      # No all combinations make sense... No need to consider spectra with a large overlap
    combinations <- combn(1:nrow(auc_overlaps), n)
    
    # Calculate and remove pre-disqualified combinations [USEFUL?]
    if(n > 6) {
      temp <- auc_overlaps; temp[lower.tri(temp)] <- 0
      dis_comb <- t(which(temp > 0.75 & temp < 2, arr.ind = TRUE))
    
      for (ii in 1:ncol(dis_comb)) {
        print(paste(ii,ncol(dis_comb),ncol(combinations)))
        combinations <- combinations[,apply(combinations, 2, function(x) !identical(dis_comb[,ii] %in% x, c(TRUE, TRUE)))]
      }
    }
    
    # Calculate total overlap for all combinations
    combi_overlap <- NULL
    combi_overlap <- apply(combinations, 2, function(x) sum(auc_overlaps[x,x]))
    
    # Find best combination
    best_combi <- combinations[,which.min(combi_overlap)]
    
    # Plot results
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
  }
  if (r == "a") {
    
  }
}
