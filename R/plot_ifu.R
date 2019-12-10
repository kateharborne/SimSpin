# Kate Harborne (last edit - 15/11/2019)
#'Plotting the flux, velocity and dipsersion image output from the IFU kinematic cube.
#'
#'The purpose of this function is to plot the useful images expracted from the IFU cube produced.
#' Input the aperture region and images produced by the \code{\link{find_lambda}} function and the
#' output will be the flux, velocity and dispersion maps.
#'
#'@param obs_data The list output from the function \code{\link{obs_data_prep}}.
#'@param obs_images The list output from the function \code{\link{obs_imgs}} containing the flux,
#' velocity and dispersion images.
#'@param reff Boolean specifying whether or not you would like the effective radius ellipse to be
#'plotted over the image. Default is FALSE.
#'@param axis_ratio The axis ratio of the effective radius ellipse. This can be taken from the
#'output of \code{\link{obs_imgs}}, or another data frame can be provided containing the semi-major
#'(\code{$a}) and semi-minor axes (\code{$b}) in kpc.
#'@param which_plots String describing which plots you wish to show - any combination of "Flux",
#'"Velocity" or "Dispersion". Default is NA, which will plot all three.
#'@return Returns three plots - a flux image, a velcoity image and a dispersion image.
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' data        = obs_data_prep(simdata = galaxy_data)
#' fluxes      = flux_grid(obs_data = data)
#' cube        = ifu_cube(obs_data  = data, flux_data = fluxes)
#' images      = obs_images(obs_data = data, ifu_datacube = cube)
#' plot_ifu(obs_data = data, obs_images = images)

plot_ifu = function(obs_data, obs_images, reff=FALSE, axis_ratio=NULL, which_plots=NA){

  ap_region      = obs_data$ap_region
  ap_region[ap_region == 0] = NA
  counts_img     = array(NA, dim=c(dim(obs_images$flux_img)[1]+20,dim(obs_images$flux_img)[2]+20))
  counts_img[11:(dim(obs_images$flux_img)[1]+10),11:(dim(obs_images$flux_img)[2]+10)] =
    (obs_images$flux_img * ap_region)
  velocity_img   = array(NA, dim=c(dim(obs_images$flux_img)[1]+20,dim(obs_images$flux_img)[2]+20))
  velocity_img[11:(dim(obs_images$flux_img)[1]+10),11:(dim(obs_images$flux_img)[2]+10)] =
    (obs_images$velocity_img * ap_region)
  dispersion_img = array(NA, dim=c(dim(obs_images$flux_img)[1]+20,dim(obs_images$flux_img)[2]+20))
  dispersion_img[11:(dim(obs_images$flux_img)[1]+10),11:(dim(obs_images$flux_img)[2]+10)] =
    obs_images$dispersion_img * ap_region

  bar_size = 5*round(((obs_data$sbinsize * obs_data$sbin) / 3)/5)

  axis_data = axis_ratio
  sbin = dim(counts_img)[1]
  sbinsize = obs_data$sbinsize
  xcen = (sbin/2)
  ycen = (sbin/2)

  if(is.na(which_plots[1])){
    which_plots = c("Flux", "Velocity", "Dispersion")
  }

  if(any(which_plots == "Flux")){
    par(family="serif", font=1, cex=1, pty="s")
    .image_nan(z = asinh(counts_img),  zlim = range(c(asinh(obs_images$flux_img))),
               col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)),
               na.color='gray', xaxt="n", yaxt="n", ann=FALSE, magmap=FALSE, family="serif", font=1)
    lines(c(5,(5+bar_size/obs_data$sbinsize)), c(5,5), col="black", lwd=2)
    points(c(5,(5+bar_size/obs_data$sbinsize)), c(5,5), col="black", lwd=3, pch="|")
    text(x = ((bar_size/obs_data$sbinsize) / 2)+5, y = 7, labels = c(paste(bar_size, " kpc")))

    fields::image.plot(legend.only = TRUE, zlim = range(c(obs_images$flux_img)),
                       col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(100)), horizontal = TRUE, family="serif", font=1,
                       legend.lab = expression("flux, Jy"))

    if (reff==TRUE){
      plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a / sbinsize, b = axis_data$b / sbinsize, angle = axis_data$ang, border="red", lwd = 5, deg=TRUE)
    }
  }

  if(any(which_plots == "Velocity")){
    par(family="serif", font=1, cex=1, pty="s")
    .image_nan(z = velocity_img,  zlim = range(c(obs_images$velocity_img)),
               col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)),
               na.color='gray', xaxt="n", yaxt="n", ann=FALSE, magmap=FALSE, family="serif", font=1)
    fields::image.plot(legend.only = TRUE, zlim = range(c(obs_images$velocity_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu"))(100)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("velocity"[LOS] * ", km s"^{-1}))
    if (reff==TRUE){
      plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a / sbinsize, b = axis_data$b / sbinsize, border="red", lwd = 5, deg=TRUE)
    }
  }

  if(any(which_plots == "Dispersion")){
    par(family="serif", font=1, cex=1, pty="s")
    .image_nan(z = dispersion_img,  zlim = range(c(obs_images$dispersion_img)),
               col=rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(200)),
               na.color='gray', xaxt="n", yaxt="n", ann=FALSE, magmap=FALSE, family="serif", font=1)
    fields::image.plot(legend.only = TRUE, zlim = range(c(obs_images$dispersion_img)), col = rev(colorRampPalette(RColorBrewer::brewer.pal(9, "RdYlBu")[1:5])(200)), horizontal = TRUE, family="serif", font=1, legend.lab = expression("dispersion"[LOS] * ", km s"^{-1}))
    if (reff==TRUE){
      plotrix::draw.ellipse(x = xcen, y = ycen, a = axis_data$a / sbinsize, b = axis_data$b / sbinsize, border="red", lwd = 5, deg=TRUE)
    }
  }
}
