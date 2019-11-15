# Kate Harborne (last edit - 23/04/2018)
#'Calculating the observed kinematics, \eqn{\lambda_R} and V/\eqn{\sigma}
#'
#'The purpose of this function is to calculate the spin parameter, \eqn{\lambda_R},
#'and ratio V/\eqn{\sigma} that would be observed given an IFU data cube. You can either
#'supply the cube created by the \code{\link{ifu_cube}} function directly, or the blurred
#'cube created by \code{\link{blur_cube}}.
#'
#'@param obs_data The list output from the function \code{\link{obs_data_prep}}.
#'@param obs_images The list output from the function \code{\link{obs_imgs}} containing the flux,
#' velocity and dispersion images.
#'@param axis_ratio The axis ratio of the effective radius ellipse. This can be taken from the
#'output of \code{\link{obs_imgs}}, or another data frame can be provided containing the semi-major
#'(\code{$a}) and semi-minor axes (\code{$b}) in kpc.
#'@return Returns a list that contains:
#' \item{\code{$obs_lambdar}}{The observed spin parameter \eqn{\lambda_R}.}
#' \item{\code{$obs_vsigma}}{The observed V/\eqn{\sigma}.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#' images      = obs_images(obs_data = data, ifu_datacube = cube)
#' kinematics  = kin_calc(obs_data = data, obs_images = images, axis_ratio = images$axis_ratio)

kin_calc = function(obs_data, obs_images, axis_ratio){

  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  vbin = obs_data$vbin # depth of cube
  sbinsize = obs_data$sbinsize

  calcregion_reff = matrix(data = 0, nrow=sbin, ncol=sbin)
  # values within reff
  radius          = matrix(data = 0, nrow=sbin, ncol=sbin)
  # radial positions within calcregion

  xcentre = sbin/2 + 0.5
  ycentre = sbin/2 + 0.5                                 # finding the centre pixel of the image
  a = axis_ratio$a / sbinsize
  b = axis_ratio$b / sbinsize

  if (a>(sbin/2)){cat("WARNING: reff > aperture, the value of $obs_lambdar produced will not be the true value evaluated at reff.", "\n")}

  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = abs(x - xcentre)
      yy = abs(y - ycentre)
      rr = ((xx^2 / a^2) + (yy^2 / b^2))
      if (rr <= 1){
        calcregion_reff[x,y] = 1
        radius[x,y] = sqrt(xx^2 + yy^2) * sbinsize
      }
    }
  }

  counts = obs_images$flux_img * calcregion_reff
  velocity = obs_images$velocity_img * calcregion_reff
  standard_dev = obs_images$dispersion_img * calcregion_reff

  lambda = sum(counts*radius*abs(velocity))/sum(counts*radius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
  vsigma = sum(counts*velocity*velocity)/sum(counts*standard_dev*standard_dev)

  output = list("obs_lr"     = lambda,
                "obs_vsigma" = sqrt(vsigma))

  return(output)

}
