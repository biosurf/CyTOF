#' Wavelength to RGB
#'
#' This function converts a given wavelength of light to an 
#' approximate RGB color value. 
#'
#' @param wavelength A wavelength value, in nanometers, in the human visual range from 380 nm through 750 nm.
#'        These correspond to frequencies in the range from 789 THz through 400 THz.
#' @param gamma The \eqn{\gamma} correction for a given display device. The linear RGB values will require
#'        gamma correction if the display device has nonlinear response.
#' @return a color string, corresponding to the result of \code{\link[grDevices]{rgb}} on the
#'        calculated R, G, B components. 
#' @source Original code taken from  http://www.noah.org/wiki/Wavelength_to_RGB_in_Python
#' @references http://stackoverflow.com/questions/1472514/convert-light-frequency-to-rgb
#' @references http://www.fourmilab.ch/documents/specrend/


wavelength_to_rgb <- function(wavelength, gamma=0.8){
  
  #
  #    Based on code by Dan Bruton
  #    http://www.physics.sfasu.edu/astro/color/spectra.html
  #    '''
  
  if (wavelength >= 380 & wavelength <= 440) {
    attenuation = 0.3 + 0.7 * (wavelength - 380) / (440 - 380)
    R = ((-(wavelength - 440) / (440 - 380)) * attenuation) ^ gamma
    G = 0.0
    B = (1.0 * attenuation) ^ gamma
  }
  else if (wavelength >= 440 & wavelength <= 490) {
    R = 0.0
    G = ((wavelength - 440) / (490 - 440)) ^ gamma
    B = 1.0
  }
  else if (wavelength >= 490 & wavelength <= 510) {
    R = 0.0
    G = 1.0
    B = (-(wavelength - 510) / (510 - 490)) ^ gamma
  }
  else if (wavelength >= 510 & wavelength <= 580) {
    R = ((wavelength - 510) / (580 - 510)) ^ gamma
    G = 1.0
    B = 0.0
  }
  else if (wavelength >= 580 & wavelength <= 645) {
    R = 1.0
    G = (-(wavelength - 645) / (645 - 580)) ^ gamma
    B = 0.0
  }
  else if (wavelength >= 645 & wavelength <= 750) {
    attenuation = 0.3 + 0.7 * (750 - wavelength) / (750 - 645)
    R = (1.0 * attenuation) ^ gamma
    G = 0.0
    B = 0.0
  }
  else {
    R = 0.0
    G = 0.0
    B = 0.0
  }
  R = R * 255
  G = G * 255
  B = B * 255
  return (rgb(floor(R), floor(G), floor(B), max=255))
}