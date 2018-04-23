# Kate Harborne (last edit - 23/04/2018)
#'Plotting the flux, velocity and dipsersion image output from the IFU kinematic cube.
#'
#'The purpose of this function is to plot the useful images expracted from the IFU cube produced.
#' Input the apperture region and images produced by the \code{\link{obs_lambda}} function and the
#' output will be the flux, velocity and dispersion maps.
#'
#'@param lambda The list output from the function \code{\link{obs_lambda}} containing the flux,
#' velocity and dispersion images.
#'@param appregion The aperture region output from the \code{\link{obs_data_prep}} function.
#'@return Returns three plots - a flux image, a velcoity image and a dispersion image.
#'@examples
#' data      = obs_data_prep(filename = system.file("extdata", 'S0_vignette', package="SimSpin"))
#' ifucube   = ifu_cube(obs_data = data)
#' reff_data = find_reff(filename     = system.file("extdata", 'S0_vignette', package="SimSpin"),
#'                       ptype        = NA,
#'                       r200         = 10,
#'                       inc_deg      = 0,
#'                       axis_ratio   = ifucube$axis_ratio,
#'                       angular_size = data$angular_size)
#' lambda_obs = obs_lambda(ifu_datacube   = ifucube,
#'                         reff_axisratio = reff_data,
#'                         sbinsize       = data$sbinsize)
#'
#' plot_ifu(lambda_obs, data$appregion)
#'

plot_ifu = function(lambda, appregion){

  appregion[appregion == 0] = NA
  counts_img     = lambda$counts_img * appregion
  velocity_img   = lambda$velocity_img * appregion
  dispersion_img = lambda$dispersion_img * appregion
  reff_ellipse   = lambda$reff_elli

  par(mfcol=c(1,1), family="serif", font=1, cex=1.1)
  magicaxis::magimage(asinh(counts_img), xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), magmap=FALSE, zlim = range(c(asinh(lambda$counts_img))), family="serif", font=1)
  fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$counts_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("flux, 10"^{-16} * "erg s"^{-1} * "cm"^{-2} * "arcsec"^{-2}))
  lines(reff_ellipse[,1], reff_ellipse[,2], col = "red", lwd = 5)

  par(family="serif", font=1, cex=1.1)
  magicaxis::magimage(velocity_img,  xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), magmap=FALSE, scale="linear")
  fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$velocity_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("velocity"[LOS] * ", km s"^{-1}))
  lines(reff_ellipse[,1], reff_ellipse[,2], col = "red", lwd = 5)

  par(family="serif", font=1, cex=1.1)
  magicaxis::magimage(dispersion_img,  xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), magmap=FALSE, scale="linear")
  fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$dispersion_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("dispersion"[LOS] * ", km s"^{-1}))
  lines(reff_ellipse[,1], reff_ellipse[,2], col = "red", lwd = 5)
}
