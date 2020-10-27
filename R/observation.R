# Author: Kate Harborne
# Date: 26/10/2020
# Title: Observation function to generate observation properties

observation = function(telescope, observing_strategy){

  ang_size      = celestial::cosdistAngScale(observing_strategy$z, ref="Planck") # angular size given z, kpc/"
  aperture_size = ang_size * telescope$fov                 # diameter size of the telescope, kpc
  sbin_size     = aperture_size / telescope$sbin           # spatial bin size (kpc per bin)
  sbin_seq      = seq(-(telescope$sbin * sbin_size) / 2,
                      (telescope$sbin * sbin_size) / 2, by=sbin_size) # spatial bin break positions
  wave_seq      = seq(telescope$wave_range[1], telescope$wave_range[2], by = telescope$wave_res) # wavelength bin break positions
  pixel_index   = seq(1,telescope$sbin*telescope$sbin, by=1)

  if (telescope$aperture_shape == "circular"){
    aperture_region = .circular_ap(telescope$sbin) * pixel_index
  }  # circular apperture mask
  if (telescope$aperture_shape == "square"){
    aperture_region = matrix(data = 1, ncol = telescope$sbin, nrow = telescope$sbin) * pixel_index
  }    # square apperture mask
  if (telescope$aperture_shape == "hexagonal"){
    aperture_region = .hexagonal_ap(telescope$sbin) * pixel_index
  }

  if (observing_strategy$blur){
    psf_fwhm = observing_strategy$fwhm
    fwhm_scaled = (psf_fwhm * ang_size)/ sbin_size  # the fwhm scaled to image pixel dimensions
    if (telescope$sbin < 25 && telescope$sbin %% 2 != 0){
      psf_dim = telescope$sbin
    } else if (telescope$sbin < 25 && telescope$sbin %% 2 == 0){
      psf_dim = telescope$sbin-1
      } else {psf_dim = 25} # dimensions of the PSF kernel

    if (observing_strategy$psf == "Gaussian"){
      psf_k = ProFit::profitMakeGaussianPSF(fwhm = fwhm_scaled, dim = c(psf_dim,psf_dim)) # the psf kernel
    }
    if (observing_strategy$psf == "Moffat"){
      psf_k = ProFit::profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = c(psf_dim,psf_dim))
    }
  } else {
    psf_fwhm = 0
    psf_k = NULL
  }

  output = list(ang_size        = ang_size,
                aperture_region = aperture_region,
                aperture_shape  = telescope$aperture_shape,
                aperture_size   = aperture_size,
                fov             = telescope$fov,
                inc_deg         = observing_strategy$inc_deg,
                inc_rad         = observing_strategy$inc_deg * (pi/180),
                lsf_fwhm        = telescope$lsf_fwhm,
                psf_fwhm        = psf_fwhm,
                psf_kernel      = psf_k,
                sbin            = telescope$sbin,
                sbin_seq        = sbin_seq,
                sbin_size       = sbin_size,
                spatial_res     = telescope$spatial_res,
                wave_bin        = length(wave_seq),
                wave_res        = telescope$wave_res,
                wave_seq        = wave_seq,
                c               = 299792.458)

  class(output) <- "observation"

  return(output)
}

.circular_ap=function(sbin){
  ap_region = matrix(data = NA, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  x = matrix(data = rep(seq(1,sbin), each=sbin), nrow = sbin, ncol = sbin)
  y = matrix(data = rep(seq(sbin,1), sbin), nrow = sbin, ncol = sbin)
  xx = x - xcentre; yy = y - ycentre
  rr = sqrt(xx^2 + yy^2)
  ap_region[rr<= sbin/2] = 1
  return(as.vector(ap_region))
}

.hexagonal_ap=function(sbin){
  ap_region = matrix(data = NA, ncol = sbin, nrow = sbin)# empty matrix for aperture mask
  xcentre = sbin/2 + 0.5; ycentre = sbin/2 + 0.5
  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = x - xcentre
      yy = y - ycentre
      rr = (2 * (sbin / 4) * (sbin * sqrt(3) / 4)) - ((sbin / 4) ) * abs(yy) - ((sbin * sqrt(3) / 4)) * abs(xx)
      if ((rr >= 0) && (abs(xx) < sbin/2) && (abs(yy) < (sbin  * sqrt(3) / 4))){
        ap_region[x,y] = 1
      }
    }
  }
  return(as.vector(ap_region))
}
