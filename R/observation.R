# Author: Kate Harborne
# Date: 26/10/2020
# Title: Observation function to generate observation properties
#
#'A class to describe the basic properties of the observation
#'
#'The purpose of this function is to generate an object that describes how the
#' galaxy has been observed.
#'
#'@param telescope A \code{telescope} object. See
#' \code{\link{telescope}} help for more details.
#'@param observing_strategy An \code{observing_strategy} object. See
#' \code{\link{observing_strategy}} help for more details.
#'@param method String to describe whether cubes output are "spectral", "gas",
#' "sf gas" or "velocity" (as in SimSpin v1) along the z-axis.
#'@return Returns an object of class "observation" that summarises the
#' properties of the observation. Run within the \code{build_datacube()}
#' function.
#'@examples
#'sami = telescope(type="SAMI")
#'conditions = observing_strategy()
#'observation_summary = observation(telescope=sami,
#'                                  observing_strategy=conditions,
#'                                  method = "spectral")
#'
observation = function(telescope, observing_strategy, method){

  ang_size      = kpc_per_arcsec(observing_strategy$distance) # angular size given z, kpc/"
  lum_dist      = Mpc(observing_strategy$distance)            # computing Luminosity Distance in units of Mpc
  aperture_size = ang_size * telescope$fov           # diameter size of the telescope, kpc
  sbin_size     = aperture_size / telescope$sbin     # spatial bin size (kpc per bin)
  sbin_seq      = seq(-(telescope$sbin * sbin_size) / 2,
                      (telescope$sbin * sbin_size) / 2, by=sbin_size) # spatial bin break positions
  wave_seq      = seq(telescope$wave_range[1], telescope$wave_range[2], by = telescope$wave_res) # wavelength bin break positions
  pixel_index   = seq(1,telescope$sbin*telescope$sbin, by=1)
  vbin_size     = min(diff(log(wave_seq))) * .speed_of_light # km/s per velocity bin
  vbin_error    = ((telescope$lsf_fwhm / telescope$wave_centre) * .speed_of_light) / (2 * sqrt(2 * log(2))) # velocity uncertainty standard deviation

  if (telescope$aperture_shape == "circular"){
    aperture_region = .circular_ap(telescope$sbin)
  }  # circular apperture mask
  if (telescope$aperture_shape == "square"){
    aperture_region = matrix(data = 1, ncol = telescope$sbin, nrow = telescope$sbin)
  }    # square apperture mask
  if (telescope$aperture_shape == "hexagonal"){
    aperture_region = .hexagonal_ap(telescope$sbin)
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
      psf_k = .gaussian_kernel(m = psf_dim, n = psf_dim, sigma = fwhm_scaled/(2*sqrt(2*log(2)))) # the psf kernel
    }
    if (observing_strategy$psf == "Moffat"){
      psf_k = ProFit::profitCubaMoffat(fwhm = fwhm_scaled, mag = 1, con = 5, dim = c(psf_dim,psf_dim))
    }
  } else {
    psf_fwhm = 0
    psf_k = NULL
  }

  wave_edges = seq(wave_seq[1] - (telescope$wave_res/2), tail(wave_seq, n=1) + (telescope$wave_res/2), by=telescope$wave_res)

  output = list(ang_size        = ang_size,
                aperture_region = aperture_region,
                aperture_shape  = telescope$aperture_shape,
                aperture_size   = aperture_size,
                date            = as.character(Sys.time()),
                fov             = telescope$fov,
                filter_name     = telescope$filter_name,
                filter          = telescope$filter,
                inc_deg         = observing_strategy$inc_deg,
                inc_rad         = observing_strategy$inc_deg * (pi/180),
                twist_deg       = observing_strategy$twist_deg,
                twist_rad       = observing_strategy$twist_deg * (pi/180),
                lsf_fwhm        = telescope$lsf_fwhm,
                lum_dist        = lum_dist,
                mass_flag       = FALSE,
                method          = method,
                origin          = paste0("SimSpin_v", packageVersion("SimSpin")),
                particle_limit  = 0,
                pointing_kpc    = xy_kpc(observing_strategy$pointing),
                pointing_deg    = xy_deg(observing_strategy$pointing),
                pixel_region    = aperture_region * pixel_index,
                psf_fwhm        = psf_fwhm,
                psf_kernel      = psf_k,
                sbin            = telescope$sbin,
                sbin_seq        = sbin_seq,
                sbin_size       = sbin_size,
                spatial_res     = telescope$spatial_res,
                signal_to_noise = telescope$signal_to_noise,
                wave_bin        = length(wave_seq),
                wave_centre     = telescope$wave_centre,
                wave_res        = telescope$wave_res,
                wave_seq        = wave_seq,
                wave_edges      = wave_edges,
                vbin_size       = vbin_size,
                vbin_error      = vbin_error,
                z               = z(observing_strategy$distance))

  return(output)
}


