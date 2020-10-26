# Kate Harborne - 16/10/2020
# SimSpin v2.0.0 - telescope class

telescope = function(type="IFU", fov=15, aperture_shape="circular", wave_range=c(3700,5700),
                     spatial_res=0.5, wave_res=1.04, lsf_fwhm=2.65){

  if(stringr::str_to_upper(type) == "SAMI"){
    output = list(type           = "SAMI",
                  fov            = 15,
                  aperture_shape = "circular",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.65,
                  sbin           = floor(15 / 0.5))
  }

  if(stringr::str_to_upper(type) == "MANGA"){
    output = list(type           = "MaNGA",
                  fov            = 22,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.5,
                  wave_res       = 1.04,
                  lsf_fwhm       = 2.8,
                  sbin           = floor(22 / 0.25))
  }

  if(stringr::str_to_upper(type) == "HECTOR"){
    output = list(type           = "Hector",
                  fov            = 30,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 0.05,
                  wave_res       = 1.6,
                  lsf_fwhm       = 1.3,
                  sbin           = floor(30 / 0.05))

  }

  if(stringr::str_to_upper(type) == "CALIFA"){
    output = list(type           = "CALIFA",
                  fov            = 30,
                  aperture_shape = "hexagonal",
                  wave_range     = c(3700,5700),
                  spatial_res    = 1,
                  wave_res       = 2,
                  lsf_fwhm       = 5.65,
                  sbin           = floor(30 / 1))

  }

  if(stringr::str_to_upper(type) == "IFU"){
    output = list(type           = "IFU",
                  fov            = fov,
                  aperture_shape = stringr::str_to_lower(aperture_shape),
                  wave_range     = wave_range,
                  spatial_res    = spatial_res,
                  wave_res       = wave_res,
                  lsf_fwhm       = lsf_fwhm,
                  sbin           = floor(fov / spatial_res))
  }

  class(output) <- "telescope"

  return(output)
}
