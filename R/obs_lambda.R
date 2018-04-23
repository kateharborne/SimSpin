# Kate Harborne (last edit - 23/04/2018)
#'Calculating the observable spin parameter, \eqn{\lambda}_R.
#'
#'The purpose of this function is to calculate the spin parameter that would be observed given an
#' IFU data cube. You can either supply the cube created by the \code{\link{ifu_cube}} function
#' directly, or the blurred cube created by \code{\link{blur_cube}}.
#'
#'@param ifu_datacube The list output from the function \code{\link{ifu_cube}} containing the mock
#' IFU cube and the apperture region image (\code{$appregion}).
#'@param reff_axisratio The semi-major and semi-minor axes output from the \code{\link{find_reff}}
#' function.
#'@param sbinsize The size of each spatial bin in kpc, output from the function
#' \code{\link{obs_data_prep}}.
#'@param dispersion_analysis (Optional) If specified as \code{TRUE}, the code will output the mean
#' and median values of the LOS velocity dispersion. Default is \code{FALSE}.
#'@return Returns a list that contains:
#' \item{\code{$obs_lambdar}}{The observed spin parameter \eqn{\lambda_R}}
#' \item{\code{$counts_img}}{The observed flux image.}
#' \item{\code{$velocity_img}}{The observed line-of-sight velocity image.}
#' \item{\code{$dispersion_img}}{The observed line-of-sight velocity dispersion image.}
#' \item{\code{$reff_ellipse}}{The coordinates of the effective radius ellipse within which
#'  \eqn{\lambda}_R is measured}
#'@examples
#'  data      = obs_data_prep(filename = system.file("extdata", 'S0_vignette', package="SimSpin"))
#'  ifucube   = ifu_cube(obs_data = data)
#'  reff_data = find_reff(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                        ptype        = NA,
#'                        r200         = 10,
#'                        inc_deg      = 0,
#'                        axis_ratio   = ifucube$axis_ratio,
#'                        angular_size = data$angular_size)
#'
#'  obs_lambda(ifu_datacube   = ifucube,
#'             reff_axisratio = reff_data,
#'             sbinsize       = data$sbinsize)
#'

obs_lambda = function(ifu_datacube, reff_axisratio, sbinsize, dispersion_analysis = FALSE){

    sbin            = length(ifu_datacube$xbin_labels)     # number of spatial bins in data cube
    vbin            = length(ifu_datacube$vbin_labels)
    calcregion_reff = array(data = rep(0,(sbin*sbin*vbin)), dim=c(sbin,sbin,vbin))
                                                           # calculating lambdaR within reff
    radius          = matrix(data = 0, nrow=sbin, ncol=sbin)
                                                           # radial positions within calcregion

    xcentre = sbin/2 + 0.5
    ycentre = sbin/2 + 0.5                                 # finding the centre pixel of the image
    a = reff_axisratio$a_kpc / sbinsize
    b = reff_axisratio$b_kpc / sbinsize

    for (x in 1:sbin){
      for (y in 1:sbin){
        xx = x - xcentre
        yy = y - ycentre
        rr = (xx^2 / a^2) + (yy^2 / b^2)
        if (rr <= 1){
          calcregion_reff[x,y,] = 1
          radius[x,y]           = sqrt(xx^2 + yy^2) * sbinsize
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

    lambda = sum(counts*radius*abs(velocity))/sum(counts*radius*(sqrt(velocity*velocity + standard_dev*standard_dev)))

    elli_x = seq(-a, a, length.out = 500)
    elli_y = (b / a) * sqrt(a^2 - elli_x^2)
    reff_elli = matrix(data=NA, nrow = 1000, ncol = 2)
    reff_elli[1:1000,1] = c(elli_x, rev(elli_x)) + rep(sbin/2, 1000)
    reff_elli[1:1000,2] = c(elli_y, rev(elli_y)*-1) + rep(sbin/2, 1000)

    if (max(reff_elli)>sbin){cat("WARNING: reff > aperture, the value of $obs_lambdar produced will not be the true value evaluated at reff.", "\n")}

    if (dispersion_analysis == TRUE) {
      mean_dispersion = mean(standard_dev)
      med_dispersion = median(standard_dev)
      dispersion_analysis = list("mean_dispersion" = mean_dispersion,
                                 "median_dispersion" = med_dispersion)
      output = list("obs_lambdar"          = lambda,
                    "counts_img"           = counts_img,
                    "velocity_img"         = velocity_img,
                    "dispersion_img"       = dispersion_img,
                    "reff_ellipse"         = reff_elli,
                    "dispersion_analysis"  = dispersion_analysis)
      } else {
        output = list("obs_lambdar"    = lambda,
                      "counts_img"     = counts_img,
                      "velocity_img"   = velocity_img,
                      "dispersion_img" = dispersion_img,
                      "reff_ellipse"   = reff_elli)
      }

  return(output)

}
