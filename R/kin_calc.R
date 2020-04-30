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
#'@param radius_type The method of computing radii - "Circular" i.e. \eqn{r^{2} = x^{2} + y{2}} or
#'"Elliptical" where r is the semi-major axis of the ellipse having an axis ratio \eqn{b/a} on
#'which the pixel lies, i.e.
#'\eqn{r^{2} = \frac{x^{2} (1 - \epsilon)^{2} + y^{2}}{(1 - \epsilon)^2}}. Default is "Both" such
#'that both \eqn{\lambda_R} values are returned.
#'@return Returns a list that contains:
#' \item{\code{$obs_lambdar}}{The observed spin parameter \eqn{\lambda_R} measured with circular
#' radii. \emph{(When \code{radius_type = "Both"} or \code{"Circular"}.)}}
#' \item{\code{$obs_elambdar}}{The observed spin parameter \eqn{\lambda_R} measured with elliptical
#' radii. \emph{(When \code{radius_type = "Both"} or \code{"Elliptical"}.)}}
#' \item{\code{$obs_vsigma}}{The observed V/\eqn{\sigma} value.}
#' \item{\code{$counts_img}}{The observed flux image.}
#' \item{\code{$velocity_img}}{The observed line-of-sight velocity image.}
#' \item{\code{$dispersion_img}}{The observed line-of-sight velocity dispersion image.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data, multi_thread=FALSE)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#' images      = obs_imgs(obs_data = data, ifu_datacube = cube)
#' kinematics  = kin_calc(obs_data = data, obs_images = images, axis_ratio = images$axis_ratio)

kin_calc = function(obs_data, obs_images, axis_ratio, radius_type = "Both"){

  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  vbin = obs_data$vbin # depth of cube
  sbinsize = obs_data$sbinsize

  calcregion_reff = matrix(data = 0, nrow=sbin, ncol=sbin)
  # values within reff
  radius          = matrix(data = 0, nrow=sbin, ncol=sbin)
  eradius         = matrix(data = 0, nrow=sbin, ncol=sbin)
  # radial positions within calcregion

  xcentre = sbin/2 + 0.5
  ycentre = sbin/2 + 0.5
  a = axis_ratio$a / sbinsize
  b = axis_ratio$b / sbinsize
  ang = axis_ratio$ang
  ang_rad = axis_ratio$ang * (pi/180)
  ellipticity = 1 - sqrt(b^2 / a^2)

  if (a*cos(ang_rad)>(sbin/2)){cat("WARNING: reff > aperture, the value of $obs_lambdar produced will not be the true value evaluated at reff.", "\n")}

  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = (x-xcentre)
      yy = (y-ycentre)
      xx = abs(x - xcentre)
      yy = abs(y - ycentre)
      rr = .ellipse(a = a, b = b, x = xx, y = yy, ang_deg = ang)
      if (rr <= 1){
        calcregion_reff[x,y] = 1
        radius[x,y] = sqrt(xx^2 + yy^2) * sbinsize
        eradius[x,y] = sqrt(((xx^2 * (1 - ellipticity)^2) + yy^2)/(1 - ellipticity)^2) * sbinsize
      }
    }
  }

  counts = obs_images$flux_img * calcregion_reff
  velocity = obs_images$velocity_img * calcregion_reff
  standard_dev = obs_images$dispersion_img * calcregion_reff

  vsigma = sum(counts*velocity*velocity)/sum(counts*standard_dev*standard_dev)

  if (radius_type == "Both" | radius_type == "both"){
    lambda = sum(counts*radius*abs(velocity))/sum(counts*radius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    elambda = sum(counts*eradius*abs(velocity))/sum(counts*eradius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    output = list("obs_lambdar"    = lambda,
                  "obs_elambdar"   = elambda,
                  "obs_vsigma"     = sqrt(vsigma),
                  "flux_img"       = counts,
                  "velocity_img"   = velocity,
                  "dispersion_img" = standard_dev)
  } else if (radius_type == "Elliptical" | radius_type == "elliptical") {
    elambda = sum(counts*eradius*abs(velocity))/sum(counts*eradius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    output = list("obs_elambdar"   = elambda,
                  "obs_vsigma"     = sqrt(vsigma),
                  "flux_img"       = counts,
                  "velocity_img"   = velocity,
                  "dispersion_img" = standard_dev)
  } else if (radius_type == "Circular" | radius_type == "circular"){
    lambda = sum(counts*radius*abs(velocity))/sum(counts*radius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    output = list("obs_lambdar"    = lambda,
                  "obs_vsigma"     = sqrt(vsigma),
                  "flux_img"       = counts,
                  "velocity_img"   = velocity,
                  "dispersion_img" = standard_dev)
  }

  return(output)

}
