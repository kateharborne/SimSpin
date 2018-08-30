# Kate Harborne (last edit - 23/04/2018)
#'Plotting the flux, velocity and dipsersion image output from the IFU kinematic cube.
#'
#'The purpose of this function is to plot the useful images expracted from the IFU cube produced.
#' Input the aperture region and images produced by the \code{\link{find_lambda}} function and the
#' output will be the flux, velocity and dispersion maps.
#'
#'@param lambda The list output from the function \code{\link{find_lambda}} containing the flux,
#' velocity and dispersion images.
#'@return Returns three plots - a flux image, a velcoity image and a dispersion image.
#'@examples
#' data      = obs_data_prep(filename = system.file("extdata", 'S0_vignette', package="SimSpin"))
#' ifucube   = ifu_cube(obs_data = data)
#' reff_data = find_props(ifu_datacube = ifucube,
#'                        sbinsize = data$sbinsize,
#'                        angular_size = data$angular_size)
#' lambda_obs = obs_lambda(ifu_datacube = ifucube,
#'                         reff_data    = reff_data,
#'                         sbinsize     = data$sbinsize)
#'
#' plot_ifu(lambda = lambda_obs, appregion = data$appregion)
#'

plot_ifu = function(lambda){

  appregion      = lambda$appregion
  appregion[appregion == 0] = NA
  counts_img     = lambda$counts_img * appregion
  velocity_img   = lambda$velocity_img * appregion
  dispersion_img = lambda$dispersion_img * appregion
  axis_data      = lambda$axis_ratio
  sbin = dim(counts_img)[1]
  sbinsize = lambda$sbinsize
  xcen = (sbin/2)
  ycen = (sbin/2)


  par(mfcol=c(1,1), family="serif", font=1, cex=1.1)
  magicaxis::magimage(asinh(counts_img), xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), magmap=FALSE, zlim = range(c(asinh(lambda$counts_img))), family="serif", font=1)
  fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$counts_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("flux, 10"^{-16} * "erg s"^{-1} * "cm"^{-2} * "arcsec"^{-2}))
  plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a_kpc / sbinsize, b = axis_data$b_kpc / sbinsize, angle = axis_data$angle - 90, border="red", lwd = 5, deg=TRUE)

  par(family="serif", font=1, cex=1.1)
  magicaxis::magimage(velocity_img,  xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), magmap=FALSE, scale="linear")
  fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$velocity_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("velocity"[LOS] * ", km s"^{-1}))
  plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a_kpc / sbinsize, b = axis_data$b_kpc / sbinsize, angle = axis_data$angle - 90, border="red", lwd = 5, deg=TRUE)

  par(family="serif", font=1, cex=1.1)
  magicaxis::magimage(dispersion_img,  xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), magmap=FALSE, scale="linear")
  fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$dispersion_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("dispersion"[LOS] * ", km s"^{-1}))
  plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a_kpc / sbinsize, b = axis_data$b_kpc / sbinsize, angle = axis_data$angle - 90, border="red", lwd = 5, deg=TRUE)
}
