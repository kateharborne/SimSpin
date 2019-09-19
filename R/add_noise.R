# Adding noise to observations
# Kate Harborne - last edit 16/09/2019

add_noise = function(ifu_datacube, noise){

  cube_wnoise = lapply(seq(1, dim(ifu_datacube$cube)[3]), function(x){ifu_datacube$cube[,,x] +
      .rtnorm(dim(ifu_datacube$cube)[1]^2, mean=ifu_datacube$cube[,,x], sd=noise, a=-ifu_datacube$cube[,,x])})
  cube_wnoise = array(unlist(cube_wnoise), dim = dim(ifu_datacube$cube)) * ifu_datacube$ap_region

  noise_cube = cube_wnoise - ifu_datacube$cube

  output = list("cube" = cube_wnoise,
                "ap_region" = ifu_datacube$ap_region,
                "xbin_labels" = ifu_datacube$xbin_labels,
                "ybin_labels" = ifu_datacube$ybin_labels,
                "vbin_labels" = ifu_datacube$vbin_labels,
                "axis_ratio"  = data.frame("a" = ifu_datacube$axis_ratio$a,
                                           "b" = ifu_datacube$axis_ratio$b),
                "noise_pcell" = noise_cube)

  return(output)
}

