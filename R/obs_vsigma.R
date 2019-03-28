# Kate Harborne (last edit - 13/09/2017)
#'Calculating the observable V/\eqn{\sigma}
#'
#'The purpose of this function is to calculate the V/\eqn{\sigma}. You can either supply the cube
#' created by the \code{\link{ifu_cube}} function directly, or the blurred cube created by
#' \code{\link{blur_cube}}.
#'
#'@param ifu_datacube The list output from the function \code{\link{ifu_cube}} containing the mock
#' IFU cube and the apperture region image (\code{$appregion}).
#'@param reff_axisratio The semi-major and semi-minor axes output from the \code{\link{find_reff}}
#' function.
#'@param sbinsize The size of each spatial bin in kpc, output from the function
#' \code{\link{obs_data_prep}}.
#'@return Returns a list that contains:
#' \item{\code{$obs_vsigma}}{The observed V/\eqn{\sigma} value.}
#' \item{\code{$counts_img}}{The observed flux image.}
#' \item{\code{$velocity_img}}{The observed line-of-sight velocity image.}
#' \item{\code{$dispersion_img}}{The observed line-of-sight velocity dispersion image.}
#'@examples
#'  galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#'  data        = obs_data_prep(simdata = galaxy_data)
#'  ifucube   = ifu_cube(obs_data = data)
#'  reff_data = find_reff(simdata      = galaxy_data,
#'                        r200         = 10,
#'                        inc_deg      = 0,
#'                        axis_ratio   = ifucube$axis_ratio,
#'                        angular_size = data$angular_size)
#'
#'  output = obs_vsigma(ifu_datacube   = ifucube,
#'                      reff_axisratio = reff_data,
#'                      sbinsize       = data$sbinsize)
#'

obs_vsigma = function(ifu_datacube, reff_axisratio, sbinsize){

  sbin            = length(ifu_datacube$xbin_labels) # number of spatial bins in data cube
  vbin            = length(ifu_datacube$vbin_labels)
  calcregion_reff = array(data = rep(0,(sbin*sbin*vbin)), dim=c(sbin,sbin,vbin))
                                                     # for calculating lambdaR within fac * reff

  xcentre = sbin/2 + 0.5
  ycentre = sbin/2 + 0.5                             # finding the centre pixel of the image
  a = reff_axisratio$a_kpc / sbinsize
  b = reff_axisratio$b_kpc / sbinsize
  ang = (reff_axisratio$angle-90) * (pi / 180)

  for (x in 1:sbin){
    for (y in 1:sbin){
      xx = x - xcentre
      yy = y - ycentre
      rr = (((xx * cos(ang)) + (yy * sin(ang)))^2 / a^2) + (((xx * sin(ang)) - (yy * cos(ang)))^2 / b^2)
      if (rr <= 1){
        calcregion_reff[x,y,] = 1
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
  }                                                        # build the velocity & dispersion maps
  velocity[(is.na(velocity))]             = 0              # mean velocity
  velocity_img[(is.na(velocity_img))]     = 0              # mean velocity image
  standard_dev[(is.na(standard_dev))]     = 0              # velocity dispersion
  dispersion_img[(is.na(dispersion_img))] = 0              # velocity dispersion image

  vsigma = sum(counts*velocity*velocity)/sum(counts*standard_dev*standard_dev)

  if (a>(sbin/2)){cat("WARNING: reff > aperture, the value of $obs_vsigma produced will not be the true value evaluated at reff.", "\n")}

  output = list("obs_vsigma"     = sqrt(vsigma),
                "counts_img"     = counts_img,
                "velocity_img"   = velocity_img,
                "dispersion_img" = dispersion_img)

  return(output)

}
