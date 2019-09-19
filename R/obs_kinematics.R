# Kate Harborne (last edit - 23/04/2018)
#'Calculating the observed kinematics, \eqn{\lambda_R} and V/\eqn{\sigma}
#'
#'The purpose of this function is to calculate the spin parameter, \eqn{\lambda_R},
#'and ratio V/\eqn{\sigma} that would be observed given an IFU data cube. You can either
#'supply the cube created by the \code{\link{ifu_cube}} function directly, or the blurred
#'cube created by \code{\link{blur_cube}}.
#'
#'@param ifu_datacube The list output from the function \code{\link{ifu_cube}} containing the mock
#' IFU cube and the apperture region image (\code{$appregion}).
#'@param reff_axisratio The semi-major and semi-minor axes output from the \code{\link{find_reff}}
#' function.
#'@param sbinsize The size of each spatial bin in kpc, output from the function
#' \code{\link{obs_data_prep}}.
#'@param radius_type The method of computing radii - "Circular" i.e. \eqn{r^{2} = x^{2} + y{2}} or
#'"Elliptical" where r is the semi-major axis of the ellipse having an axis ratio \eqn{b/a} on
#'which the pixel lies, i.e.
#'\eqn{r^{2} = \frac{x^{2} (1 - \epsilon)^{2} + y^{2}}{(1 - \epsilon)^2}. Default is "Both" such
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
#'  galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#'  data        = obs_data_prep(simdata = galaxy_data)
#'  ifucube     = ifu_cube(obs_data = data, threshold = 20)
#'  reff_data = find_reff(simdata      = galaxy_data,
#'                        r200         = 200,
#'                        inc_deg      = 0,
#'                        axis_ratio   = ifucube$axis_ratio,
#'                        angular_size = data$angular_size)
#'
#'  output = obs_kinematics(ifu_datacube   = ifucube,
#'                          reff_axisratio = reff_data,
#'                          sbinsize       = data$sbinsize)
#'

obs_kinematics = function(ifu_datacube, reff_axisratio, sbinsize, radius_type = "Both"){

  sbin            = length(ifu_datacube$xbin_labels)       # number of spatial bins in data cube
  vbin            = length(ifu_datacube$vbin_labels)
  calcregion_reff = array(data = rep(0,(sbin*sbin*vbin)),
                          dim=c(sbin,sbin,vbin))           # calculating lambdaR within reff
  radius          = matrix(data = 0, nrow=sbin, ncol=sbin) # radial positions within calcregion
  eradius         = matrix(data = 0, nrow=sbin, ncol=sbin)

  xcentre = sbin/2 + 0.5
  ycentre = sbin/2 + 0.5                                 # finding the centre pixel of the image
  a = reff_axisratio$a_kpc / sbinsize
  b = reff_axisratio$b_kpc / sbinsize
  ang = (reff_axisratio$angle-90) * (pi / 180)
  ellipticity = 1 - sqrt(b^2 / a^2)

  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = abs(x - xcentre)
      yy = abs(y - ycentre)
      rr = (((xx * cos(ang)) + (yy * sin(ang)))^2 / a^2) + (((xx * sin(ang)) - (yy * cos(ang)))^2 / b^2)
      if (rr <= 1){
        calcregion_reff[x,y,] = 1
        radius[x,y] = sqrt(xx^2 + yy^2) * sbinsize
        eradius[x,y] = sqrt(((xx^2 * (1 - ellipticity)^2) + yy^2)/(1 - ellipticity)^2) * sbinsize
      }
    }
  }
  # creating two arrays - one of a multiplier to consider lambdaR within Reff,
  #  and one of the radial values within Reff

  cube_reff  = ifu_datacube$cube * calcregion_reff
  counts     = apply(cube_reff, c(1,2), sum)
  counts_img = apply(ifu_datacube$cube, c(1,2), sum)
  velocity       = matrix(data=0, nrow=sbin, ncol=sbin)
  velocity_img   = matrix(data=0, nrow=sbin, ncol=sbin)
  standard_dev   = matrix(data=0, nrow=sbin, ncol=sbin)
  dispersion_img = matrix(data=0, nrow=sbin, ncol=sbin)
  for (c in 1:sbin){
    for (d in 1:sbin){
      velocity[c,d]       = .meanwt(ifu_datacube$vbin_labels, cube_reff[c,d,])
      standard_dev[c,d]   = sqrt(.varwt(ifu_datacube$vbin_labels, cube_reff[c,d,], velocity[c,d]))
      velocity_img[c,d]   = .meanwt(ifu_datacube$vbin_labels, ifu_datacube$cube[c,d,])
      dispersion_img[c,d] = sqrt(.varwt(ifu_datacube$vbin_labels, ifu_datacube$cube[c,d,], velocity_img[c,d]))
    }
  }
  # building the velocity and dispersion images
  velocity[(is.na(velocity))]             = 0            # mean velocity
  velocity_img[(is.na(velocity_img))]     = 0            # mean velocity image
  standard_dev[(is.na(standard_dev))]     = 0            # velocity dispersion
  dispersion_img[(is.na(dispersion_img))] = 0            # velocity dispersion image

  vsigma = sum(counts*velocity*velocity)/sum(counts*standard_dev*standard_dev)

  if (a>(sbin/2)){cat("WARNING: reff > aperture, the value of $obs_lambdar produced will not be the true value evaluated at reff.", "\n")}
  if (radius_type == "Both" | radius_type == "both"){
    lambda = sum(counts*radius*abs(velocity))/sum(counts*radius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    elambda = sum(counts*eradius*abs(velocity))/sum(counts*eradius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    output = list("obs_lambdar"    = lambda,
                  "obs_elambdar"   = elambda,
                  "obs_vsigma"     = sqrt(vsigma),
                  "counts_img"     = counts_img,
                  "velocity_img"   = velocity_img,
                  "dispersion_img" = dispersion_img)
  } else if (radius_type == "Elliptical" | radius_type == "elliptical") {
    elambda = sum(counts*eradius*abs(velocity))/sum(counts*eradius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    output = list("obs_elambdar"   = elambda,
                  "obs_vsigma"     = sqrt(vsigma),
                  "counts_img"     = counts_img,
                  "velocity_img"   = velocity_img,
                  "dispersion_img" = dispersion_img)
  } else if (radius_type == "Circular" | radius_type == "circular"){
    lambda = sum(counts*radius*abs(velocity))/sum(counts*radius*(sqrt(velocity*velocity + standard_dev*standard_dev)))
    output = list("obs_lambdar"    = lambda,
                  "obs_vsigma"     = sqrt(vsigma),
                  "counts_img"     = counts_img,
                  "velocity_img"   = velocity_img,
                  "dispersion_img" = dispersion_img)
  }

  return(output)

}
