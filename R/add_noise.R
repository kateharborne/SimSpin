# Adding noise to observations
# Kate Harborne - last edit 16/09/2019

add_noise = function(image_grid, SN=50, mag_threshold=25.2, mag_zero=0, pixel_sscale=0.2){
  signal = sum(image_grid$flux_img)
  s0 = rnorm(dim(image_grid$flux_img)[1]^2, sd=sqrt(signal))
  skyRMS = ProFound::profoundSB2Flux(mag_threshold, mag_zero, pixel_sscale)
  ss = rnorm(dim(image_grid$flux_img)[1]^2, sd=skyRMS)

  if (signal/sum(s0 + ss) > SN){
    sd = rpois(dim(image_grid$flux_img)[1]^2, as.vector(image_grid$flux_img))
  }


  counts_img = image_grid$flux_img + s0 + ss

}

