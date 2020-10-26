# Kate Harborne - 16/10/2020
# SimSpin v2.0.0 - observing_strategy class

observing_strategy = function(z = 0.1, inc_deg = 70, blur = F, fwhm=2, psf="Gaussian"){

  if(blur){
    if(stringr::str_to_upper(psf) == "GAUSSIAN"){
      output = list(z       = z,
                    inc_deg = inc_deg,
                    blur    = T,
                    fwhm    = fwhm,
                    psf     = "Gaussian")
    }
    if(stringr::str_to_upper(psf) == "MOFFAT"){
      output = list(z       = z,
                    inc_deg = inc_deg,
                    blur    = T,
                    fwhm    = fwhm,
                    psf     = "Moffat")
    }
  } else {
    output = list(z       = z,
                  inc_deg = inc_deg,
                  blur    = F)
  }

  class(output) <- "observing_strategy"

  return(output)

}
