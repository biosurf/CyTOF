# l = c(355,405,488,561,640)
# m = "all"
# c = 0.2

emission_recalc <- function(l, m, c) {
  load("~/Dropbox/Git_repos/biosurf/flow/fluorochromes_all.Rdata")
  f_sub <- list()
  peak_wl <- c()
  
  if(class(m) == "list") {
    fluorochromes_all <- fluorochromes_all[which(names(fluorochromes_all) %in% unname(unlist(m)))]
  }
  for(i in 1:length(fluorochromes_all)) {
    x <- fluorochromes_all[[i]]
    x$Emission[x$Emission < 0] <- 0
    l_max <- c()
    for(ii in 1:length(l)) {
      l_max <- append(l_max, x$Excitation[which(x$Wavelength==l[ii])])
    }
    l_max <- which.max(l_max)
    if ((x$Excitation[which(x$Wavelength==l[l_max])]/100) >= c) {
      f_sub[[names(fluorochromes_all)[i]]] <- x
      f_sub[[names(fluorochromes_all)[i]]]$Emission <- (f_sub[[names(fluorochromes_all)[i]]]$Emission * x$Excitation[which(x$Wavelength==l[l_max])]) / 100
      peak_wl <- append(peak_wl,x$Wavelength[which.max(x$Emission)])
    }
    else {print(paste(names(fluorochromes_all)[i], " omitted. Relative emission intensity is below ", c*100, "%", sep=""))}
  }
  
  f_sub <- f_sub[order(peak_wl)]
  
  return(f_sub)
}