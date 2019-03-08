# Kate Harborne (last edit - 23/04/2018)
#'Plotting the flux, velocity and dipsersion image output from the IFU kinematic cube.
#'
#'The purpose of this function is to plot the useful images expracted from the IFU cube produced.
#' Input the aperture region and images produced by the \code{\link{find_lambda}} function and the
#' output will be the flux, velocity and dispersion maps.
#'
#'@param lambda The list output from the function \code{\link{find_lambda}} containing the flux,
#' velocity and dispersion images.
#'@param reff Boolean specifying whether or not you would like the effective radius ellipse to be
#'plotted over the image. Default is TRUE.
#'@param which_plots String describing which plots you wish to show - any combination of "Flux",
#'"Velocity" or "Dispersion". Default is NA, which will plot all three.
#'@return Returns three plots - a flux image, a velcoity image and a dispersion image.
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' lambdar = find_lambda(simdata = galaxy_data)
#' plot_ifu(lambdar)
#'

plot_ifu = function(lambda, reff=TRUE, which_plots=NA){

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

  if(is.na(which_plots[1])){
    which_plots = c("Flux", "Velocity", "Dispersion")
  }

  if(any(which_plots == "Flux")){
    par(mfcol=c(1,1), family="serif", font=1, cex=1.1)
    magicaxis::magimage(asinh(counts_img), xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), magmap=FALSE, zlim = range(c(asinh(lambda$counts_img))), family="serif", font=1)
    fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$counts_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("flux, 10"^{-16} * "erg s"^{-1} * "cm"^{-2} * "arcsec"^{-2}))
    if (reff==TRUE){
      plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a_kpc / sbinsize, b = axis_data$b_kpc / sbinsize, angle = axis_data$angle - 90, border="red", lwd = 5, deg=TRUE)
    }
  }

  if(any(which_plots == "Velocity")){
    par(family="serif", font=1, cex=1.1)
    magicaxis::magimage(velocity_img,  xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), magmap=FALSE, scale="linear")
    fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$velocity_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("velocity"[LOS] * ", km s"^{-1}))
    if (reff==TRUE){
      plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a_kpc / sbinsize, b = axis_data$b_kpc / sbinsize, angle = axis_data$angle - 90, border="red", lwd = 5, deg=TRUE)
    }
  }

  if(any(which_plots == "Dispersion")){
    par(family="serif", font=1, cex=1.1)
    magicaxis::magimage(dispersion_img,  xaxt="n", yaxt="n", ann=FALSE, col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(200)), magmap=FALSE, scale="linear")
    fields::image.plot(legend.only = TRUE, zlim = range(c(lambda$dispersion_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(200)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("dispersion"[LOS] * ", km s"^{-1}))
    if (reff==TRUE){
      plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a_kpc / sbinsize, b = axis_data$b_kpc / sbinsize, angle = axis_data$angle - 90, border="red", lwd = 5, deg=TRUE)
    }
  }
}
