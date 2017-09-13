# Kate Harborne (last edit - 12/09/2017)
#'Creating a set of mock IFU observational images.
#'
#'The purpose of this function is to construct an IFU images. It accepts output parameters from \code{obs_data_prep()} and returns
#' a counts image, a velocity image and a dispersion image along with an apperture region image that describes the shape of the apperture.
#'
#'@param obs_data The list output from the \code{obs_data_prep()} function.
#'@param threshold The minimum number of counts in the image.
#'@param type \emph{Optional} The type of data returned. The option \code{type = "physical"} returns distances in units of kpc,
#' but specifying \code{type = "observer"} will return distances in arcseconds.
#'@return Returns a list containing a mock IFU data images including a counts image (\code{$counts_img}), a velocity image (\code{$velocity_img})
#' and a dispersion image (\code{$dispersion_img}) as required for calculating the observed spin parameter. Also contained is an image that
#' describes the shape of the apperture (\code{$appregion}) such that any further convolutions that are applied to mimic beam smearing or seeing
#' can be trimmed to the appropriate apperture shape.
#'@examples
#' \dontrun{
#' data = obs_data_prep()
#'
#' ifu_img(obs_data  = data,
#'         threshold = 0)
#'
#' ifu_img(obs_data  = data,
#'         threshold = 20,
#'         type      = "observer")
#' }
#'

ifu_img = function(obs_data, threshold, type="physical") {

  galaxy_obs      = obs_data$galaxy_obs # data to populate the aperture
  sbin            = obs_data$sbin # number of spatial bins along the x and z axes
  sbin_breaks     = seq(-(sbin * obs_data$sbinsize) / 2, (sbin * obs_data$sbinsize) / 2, by=obs_data$sbinsize)
  galaxy_obs$binx = cut(galaxy_obs$x, breaks=sbin_breaks, labels=F)
  galaxy_obs$binz = cut(galaxy_obs$z_obs, breaks=sbin_breaks, labels=F)
  galaxy_obs$binn = galaxy_obs$binx + (sbin * galaxy_obs$binz) - sbin # labelling data based on it's position in the aperture

  if (type == "observer"){
    xbin_labels =  xbin_labels / obs_data$angular_size
    zbin_labels =  zbin_labels / obs_data$angular_size
  }

  counts      = with(galaxy_obs, as.integer(binn))
  counts_df   = as.data.frame(table(counts))
  counts_flat = data.frame("counts"=rep(0,sbin*sbin))
  counts_flat[as.integer(as.vector(counts_df$counts)),] = as.integer(as.vector(counts_df$Freq))
  counts_img   = (array(t(as.matrix(counts_flat, ncol=1)), dim=c(sbin,sbin)))
  appregion    = counts_img
  appregion[(appregion > 0)] = 1
  counts_img[(counts_img < threshold)] = 0

  V_counts = as.data.frame(xtabs(formula = vy_obs ~ binn, data = galaxy_obs))
  V_flat   = data.frame("velocity"=rep(0,sbin*sbin))
  V_flat[as.integer(as.vector(V_counts$binn)),] = as.integer(as.vector(V_counts$Freq))
  velocity_img = (array(t(as.matrix(V_flat, ncol=1)), dim=c(sbin,sbin))) / counts_img
  velocity_img[is.nan(velocity_img) | is.infinite(velocity_img)] = 0

  V2_counts = as.data.frame(xtabs(formula = (vy_obs * vy_obs) ~ binn, data = galaxy_obs))
  V2_flat   = data.frame("velocity2"=rep(0,sbin*sbin))
  V2_flat[as.integer(as.vector(V2_counts$binn)),] = as.integer(as.vector(V2_counts$Freq))
  velocity2_img = (array(t(as.matrix(V2_flat, ncol=1)), dim=c(sbin,sbin))) / counts_img
  velocity2_img[is.nan(velocity2_img) | is.infinite(velocity2_img)] = 0
  dispersion_img = sqrt(velocity2_img - (velocity_img)^2)

  output = list("counts_img" = counts_img, "velocity_img" = velocity_img, "dispersion_img" = dispersion_img, "appregion" = appregion)

  return(output)

}
