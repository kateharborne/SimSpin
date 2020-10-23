# Kate Harborne - 16/10/2020
# SimSpin v2.0.0 - telescope class

telescope = function(type="IFU", fov=15, aperture_shape="circular", wave_range=c(3700,5700),
                     spatial_res=0.5, wave_res=1.04, lsf_fwhm=2.65){

  if(type == "SAMI" | type == "sami" | type == "Sami"){
    output = list(type           = "SAMI",
                  fov            = 15,
                  aperture_shape = "circular",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.65)
  }

  if(type == "MaNGA" | type == "manga" | type == "MANGA" | type == "Manga"){
    output = list(type           = "MaNGA",
                  fov            = 22,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.65)
  }

  if(type == "Hector" | type == "HECTOR" | type == "hector"){
    output = list(type           = "MaNGA",
                  fov            = 22,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.65)

  }

  if(type == "CALIFA" | type == "Califa" | type == "califa"){
    output = list(type           = "MaNGA",
                  fov            = 22,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.65)

  }

  if(type == "IFU" | type == "ifu" | type == "Ifu"){
    output = list(type           = "IFU",
                  fov            = fov,
                  aperture_shape = aperture_shape,
                  wave_range     = wave_range,
                  spatial_res    = spatial_res,
                  wave_res       = wave_res,
                  lsf_fwhm       = lsf_fwhm)
  }

  class(output) <- "telescope"

  return(output)
}
