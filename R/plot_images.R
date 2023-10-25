# Author: Kate Harborne
# Date: 25/10/2023
# Title: plot_images - a suite of function for plotting pretty images

#' Plotting pretty flux images
#'
#' A function to produce a plot of the flux image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param flux_image Numeric array containing the flux image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_flux(cube$observed_images$flux_image)

plot_flux <- function(flux_image, fig = c(0,1,0,1), new=F,
                      units = expression("Flux, CGS"), main="",
                      radii = NA, na.color = "white", zlim = NA, legend = T,
                      titleshift = -4, labN=5, radii_col="red", ...){

  Flux = flux_image
  im_dim = dim(Flux)/2
  flux_val = c(min(Flux, na.rm = T), max(Flux, na.rm = T))

  if (all(flux_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  flux_map_cols = cmocean::cmocean("thermal", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = Flux, zlim =  if(is.na(zlim[1])){flux_val}else{zlim}, col = flux_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){flux_val}else{zlim}, scale = c(1, 1/20),
              col = flux_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty mass images
#'
#' A function to produce a plot of the mass image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param mass_image Numeric array containing the mass image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity", mass_flag=TRUE)
#' plot_mass(cube$observed_images$mass_image)

plot_mass <- function(mass_image, fig = c(0,1,0,1), new=F,
                      units = expression("Mass, M"["sol"]), main="",
                      radii = NA, na.color = "white", zlim = NA, legend = T,
                      titleshift = -4, labN=5, radii_col="red", ...){

  Mass = mass_image
  im_dim = dim(Mass)/2
  mass_val = c(min(Mass, na.rm = T), max(Mass, na.rm = T))

  if (all(mass_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  mass_map_cols = cmocean::cmocean("thermal", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = Mass, zlim =  if(is.na(zlim[1])){mass_val}else{zlim}, col = mass_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){mass_val}else{zlim}, scale = c(1, 1/20),
               col = mass_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
               titleshift = titleshift)
  }

}



#' Plotting pretty velocity images
#'
#' A function to produce a plot of the velocity image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param velocity_image Numeric array containing the velocity image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_velocity(cube$observed_images$velocity_image)

plot_velocity <- function(velocity_image, fig = c(0,1,0,1), new=F,
                          units = expression("velocity"[LOS] * ", km s"^{-1}), main="",
                          radii = NA, na.color = "white", zlim = NA, legend = T,
                          titleshift = -4, labN=5, radii_col="red", ...){

  V = velocity_image
  im_dim = dim(V)/2
  vel_val = max(c(abs(min(V, na.rm = T)), abs(max(V, na.rm = T))))

  if (all(vel_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  velo_map_cols = cmocean::cmocean("balance", version = "2.0", start = .1, end=.9)(100)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = V, zlim = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, col = velo_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, scale = c(1, 1/20),
              col = velo_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty dispersion images
#'
#' A function to produce a plot of the dispersion image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param dispersion_image Numeric array containing the dispersion image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_dispersion(cube$observed_images$dispersion_image)

plot_dispersion <- function(dispersion_image, fig = c(0,1,0,1), new=F,
                            units = expression("dispersion"[LOS] * ", km s"^{-1}), main="",
                            radii = NA, na.color = "white", zlim = NA, legend=T,
                            titleshift = -4, labN=5, radii_col="red", ...){

  disp_map = dispersion_image
  im_dim = dim(disp_map)/2
  disp_val = c(floor(min(disp_map, na.rm=T)), max(disp_map, na.rm=T))

  if (all(disp_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  disp_map_cols = cmocean::cmocean("solar", version = "2.0", start = .1, end=.9)(100)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = disp_map, zlim = if(is.na(zlim[1])){disp_val}else{zlim},
             col = disp_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){disp_val}else{zlim}, scale = c(1, 1/20),
              col = disp_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty higher-order kinematic (h3) images
#'
#' A function to produce a plot of the h3 image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param h3_image Numeric array containing the h3 image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_h3(cube$observed_images$h3_image)

plot_h3   <- function(h3_image, fig = c(0,1,0,1), new=F,
                      units = expression("h"[3]), main="",
                      radii = NA, na.color = "white", zlim = NA, legend = T,
                      titleshift = -4, labN=5, radii_col="red", ...){

  V = h3_image
  im_dim = dim(V)/2
  vel_val = max(c(abs(min(V, na.rm = T)), abs(max(V, na.rm = T))))

  if (all(vel_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  velo_map_cols = cmocean::cmocean("balance", version = "2.0", start = .1, end=.9)(100)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = V, zlim = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, col = velo_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, scale = c(1, 1/20),
              col = velo_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty higher-order kinematic (h4) images
#'
#' A function to produce a plot of the h4 image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param h4_image Numeric array containing the h4 image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_h4(cube$observed_images$h4_image)

plot_h4   <- function(h4_image, fig = c(0,1,0,1), new=F,
                      units = expression("h"[4]), main="",
                      radii = NA, na.color = "white", zlim = NA, legend = T,
                      titleshift = -4, labN=5, radii_col="red", ...){

  V = h4_image
  im_dim = dim(V)/2
  vel_val = max(c(abs(min(V, na.rm = T)), abs(max(V, na.rm = T))))

  if (all(vel_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  velo_map_cols = cmocean::cmocean("balance", version = "2.0", start = .1, end=.9)(100)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = V, zlim = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, col = velo_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main = main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){c(-vel_val,vel_val)}else{zlim}, scale = c(1, 1/20),
              col = velo_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty age images
#'
#' A function to produce a plot of the stellar age image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param age_image Numeric array containing the age image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_age(cube$raw_images$age_image)


plot_age <- function(age_image, fig = c(0,1,0,1), new=F,
                     units = expression("Age, Gyr"), main="", radii = NA,
                     na.color = "white", zlim = NA, legend=T,
                     titleshift = -4, labN=5, radii_col="red", ...){

  age_map = age_image
  im_dim = dim(age_map)/2
  age_val = c(floor(min(age_map, na.rm=T)/2), max(age_map, na.rm=T))

  if (all(age_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  age_map_cols = cmocean::cmocean("deep", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = age_map, zlim = if(is.na(zlim[1])){age_val}else{zlim}, col = age_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){age_val}else{zlim}, scale = c(1, 1/20),
              col = age_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty metallicity images
#'
#' A function to produce a plot of the stellar metallicity image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param metallicity_image Numeric array containing the metallicity image.
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_metallicity(cube$raw_images$metallicity_image)

plot_metallicity <- function(metallicity_image, fig = c(0,1,0,1), new=F,
                             units = expression("log10(Z/Z"[solar]*")"), main="",
                             na.color = "white", zlim = NA, legend=T, radii = NA,
                             titleshift = -4, labN=5, radii_col="red", ...){

  met_map = metallicity_image
  im_dim = dim(met_map)/2
  met_val = c(floor(min(met_map, na.rm=T)/2), max(met_map, na.rm=T))

  if (all(met_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  met_map_cols = cmocean::cmocean("dense", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = met_map, zlim = if(is.na(zlim[1])){met_val}else{zlim},
             col = met_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){met_val}else{zlim}, scale = c(1, 1/20),
              col = met_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
              titleshift = titleshift)
  }

}

#' Plotting pretty particle number images
#'
#' A function to produce a plot of the particle per pixel image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param particle_image Numeric array containing the particle per pixel image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_gadget = system.file("extdata", "SimSpin_example_Gadget",
#' package = "SimSpin")
#' ss_gadget = make_simspin_file(ss_pd_gadget, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_gadget,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity")
#' plot_particles(cube$raw_images$particle_image)

plot_particles <- function(particle_image, fig = c(0,1,0,1), new=F,
                           units = expression("Number of particles"), main="", radii = NA,
                           na.color = "white", zlim = NA, legend=T,
                           titleshift = -4, labN=5, radii_col="red", ...){

  part_map = particle_image
  im_dim = dim(part_map)/2
  part_val = c(floor(min(part_map, na.rm=T)/2), max(part_map, na.rm=T))

  if (all(part_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  part_map_cols = cmocean::cmocean("deep", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = part_map, zlim = if(is.na(zlim[1])){part_val}else{zlim}, col = part_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){part_val}else{zlim}, scale = c(1, 1/20),
               col = part_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
               titleshift = titleshift)
  }

}

#' Plotting pretty gas SFR images
#'
#' A function to produce a plot of the gas star formation rate image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param SFR_image Numeric array containing the particle per pixel image
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5",
#' package = "SimSpin")
#' ss_eagle = make_simspin_file(ss_pd_eagle, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_eagle,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "gas")
#' plot_SFR(cube$raw_images$SFR_image)

plot_SFR <- function(SFR_image, fig = c(0,1,0,1), new=F,
                     units = expression("SFR, M"["sol"]*"/yr"), main="", radii = NA,
                     na.color = "white", zlim = NA, legend=T,
                     titleshift = -4, labN=5, radii_col="red",...){

  part_map = SFR_image
  im_dim = dim(part_map)/2
  part_val = c(floor(min(part_map, na.rm=T)/2), max(part_map, na.rm=T))

  if (all(part_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  part_map_cols = cmocean::cmocean("deep", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = part_map, zlim = if(is.na(zlim[1])){part_val}else{zlim}, col = part_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){part_val}else{zlim}, scale = c(1, 1/20),
               col = part_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
               titleshift = titleshift)
  }

}

#' Plotting pretty log10(O/H) images
#'
#' A function to produce a plot of the gas log10(O/H) image produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param OH_image Numeric array containing the metallicity image.
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param units String describing the units of the values contained
#' in the image.
#' @param main Image title, default "".
#' @param radii list - containing a, b, and ang (if wishing to plot half-mass
#' radii ellipse).
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param legend Boolean to determine if the colour bar axis should be printed
#' at the bottom of the image. Default is T.
#' @param titleshift Numeric. Describes the distance between the colour bar and
#' the units. Default is -4.
#' @param labN Numeric. Describes the minimum number of numeric labels added to
#' the colour bar. Default is 5.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5",
#' package = "SimSpin")
#' ss_eagle = make_simspin_file(ss_pd_eagle, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_eagle,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "gas")
#' plot_OH(cube$raw_images$OH_image)

plot_OH <- function(OH_image, fig = c(0,1,0,1), new=F,
                    units = expression("log10(O/H) + 12"), main="",
                    na.color = "white", zlim = NA, legend=T, radii = NA,
                    titleshift = -4, labN=5, radii_col="red", ...){

  met_map = OH_image
  im_dim = dim(met_map)/2
  met_val = c(floor(min(met_map, na.rm=T)/2), max(met_map, na.rm=T))

  if (all(met_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  met_map_cols = cmocean::cmocean("dense", version = "2.0", start = .1, end=1)(50)

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = met_map, zlim = if(is.na(zlim[1])){met_val}else{zlim},
             col = met_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
  if (!is.na(radii[1])){
    plotrix::draw.ellipse(im_dim[1], im_dim[2], radii$a, radii$b, radii$ang, border = radii_col, density = NULL, lwd=2)
  }
  if (legend){
    .magcolbar(position = "bottom", range = if(is.na(zlim[1])){met_val}else{zlim}, scale = c(1, 1/20),
               col = met_map_cols, orient = "h", inset = -1/20, labN=labN, title = units,
               titleshift = titleshift)
  }

}

#' Plotting pretty Voronoi bin images
#'
#' A function to produce a plot of the Voronoi bins produced by
#' \code{build_datacube()} with associated colour bar and labels.
#'
#' @param voronoi_bins Numeric array containing the metallicity image.
#' @param fig Numeric array of length 4 describing the boundary of the image
#' @param new Boolean. Should the image be added to the existing plot? Default
#' is FALSE.
#' @param main Image title, default "".
#' @param na.color String. Colour given to NA values in the image.
#' @param zlim Numeric array of length 2. Describing the numeric range of
#' colours in the image. Default is NA, in which the range will be described by
#' the minimum and maximum values in the image.
#' @param ... Further variables passed to magimage. See
#' \code{\link[magicaxis]{magimage}} for further details.
#' @return Returns an image to the plotting window of the input
#' \code{build_datacube} image.
#' @examples
#' ss_pd_eagle = system.file("extdata", "SimSpin_example_EAGLE.hdf5",
#' package = "SimSpin")
#' ss_eagle = make_simspin_file(ss_pd_eagle, write_to_file = FALSE)
#' cube = build_datacube(simspin_file = ss_eagle,
#'                       telescope = telescope(type="SAMI"),
#'                       observing_strategy = observing_strategy(),
#'                       method = "velocity", voronoi_bin = TRUE)
#' plot_voronoi(cube$raw_images$voronoi_bins)

plot_voronoi <- function(voronoi_bins, fig = c(0,1,0,1), new=F,
                         main="", na.color = "white", zlim = NA,
                         ...){

  bin_map = voronoi_bins
  im_dim = dim(bin_map)/2
  bin_val = c(min(bin_map, na.rm=T), max(bin_map, na.rm=T))

  if (all(bin_val == 0)){
    stop("Image contains only '0'. No image can be produced in this case. \n
         Please check your build_datacube function and try again.")
  }

  bin_map_cols = cmocean::cmocean("phase", version = "2.0", start = .1, end=1)(max(bin_map, na.rm=T))

  par(pty="s", fig=fig, xpd=FALSE, ps=12, cex=1, new=new); options(scipen = 1)
  .image_nan(z = bin_map, zlim = if(is.na(zlim[1])){bin_val}else{zlim},
             col = bin_map_cols, na.color = na.color, xaxt="n",
             yaxt="n", ann=FALSE, magmap=FALSE, family="mono", font=1, main=main, ...)
}

.image_nan <- function(z, zlim, col, na.color='gray', ...){
  # Function for plotting an image with NAs as a set colour
  # z        :==: array corresponding to the image to be plotted
  # zlim     :==: colour range associated with the image for
  #                the colour palette
  # col      :==: the colour palette to be used for the image
  # na.color :==: string describing the color used for NA

  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.na <- zlim[2] + zstep # new z for NA

  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na

  zlim[2] <- zlim[2] + zstep # extend top limit to include the new value na

  col <- c(col, na.color) # we construct the new color range by including: na.color

  magicaxis::magimage(z=z, zlim=zlim, col=col, ...) # we finally call image(...)
}

.magcolbar <- # function modified from magicaxis to generate colour bars
  function(position='topright',range=c(0,1),orient='v',log=FALSE,col=hcl.colors(21),scale=c(1/4,1/20), inset=1/40,labN=5,title='',titleshift=0,centrealign='rb',clip='',cex=1,...){
    usercoord=par()$usr
    xlogcheck=FALSE;ylogcheck=FALSE

    xlo=usercoord[1];xhi=usercoord[2];ylo=usercoord[3];yhi=usercoord[4]
    xdiff=xhi-xlo;ydiff=yhi-ylo

    if(orient=='h'){
      xl=xlo+xdiff/2-xdiff*scale[1]/2
      yb=ylo+ydiff/2-ydiff*scale[2]/2
      xr=xlo+xdiff/2+xdiff*scale[1]/2
      yt=ylo+ydiff/2+ydiff*scale[2]/2
      align=centrealign
      if(length(grep('bottom',position))>0){aligntext='lt'; yb=ylo+ydiff*inset;yt=ylo+ydiff*inset+ydiff*scale[2]}
    }

    rangetemp=range
    legend=magicaxis::maglab(rangetemp,labN,log=log,trim=F)

    roughNscale=(max(legend$labat)-min(legend$labat))/(rangetemp[2]-rangetemp[1])

    colremap=magicaxis::magmap(data=seq(min(legend$labat),max(legend$labat),length=length(col)*roughNscale),locut=rangetemp[1],hicut=rangetemp[2],type='num',range=c(1,length(col)),clip=clip)$map
    col=col[round(colremap,digits=0)]

    if(orient=='h'){plotrix::color.legend(xl,yb,xr,yt,legend=legend$exp,rect.col=col,cex=cex,align=align,gradient='x',...)}
    if(orient=='h' & aligntext=='lt'){text((xl+xr)/2,yt+(titleshift)*ydiff/20, labels=title,adj=c(0.5,0.5),srt=0,cex=cex, xpd=NA)}

    par(xlog=xlogcheck)
    par(ylog=ylogcheck)
    par(usr=usercoord)
  }
