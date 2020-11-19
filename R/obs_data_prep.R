# Kate Harborne (last edit - 14/11/19)
#'Prepare simulation data for observational kinematic analysis.
#'
#'The purpose of this function is to calculate the properties necessary for constructing an IFU
#' observation data cube of a simulated galaxy at a user specified inclination and redshift.
#'
#'@param simdata The simulation information data.frame output by \code{\link{sim_data}}.
#'@param r200 The virial radius specified in the simulation, kpc.
#'@param z The galaxy redshift.
#'@param fov The field of view of the IFU, diameter in arcseconds.
#'@param ap_shape The shape of the field of view, with options "circular", "square" or "hexagonal".
#'@param central_wvl The central filter wavelength used for the observation, given in angstroms.
#'@param lsf_fwhm The line spread function full-width half-max, given in angstroms.
#'@param pixel_sscale The corresponding spatial pixel scale associated with a given telescope
#' output in arcseconds.
#'@param pixel_vscale The corresponding velocity pixel scale associated with a given telescope
#' filter output in angstroms.
#'@param inc_deg The inclination at which to observe the galaxy in degrees.
#'@param align Boolean indicating whether or not to align the semi-major axis with the x-axis.
#'@return Returned is a list that contains:
#' \item{\code{$galaxy_obs}}{A data frame of the observed particle information (which contains the
#'   galaxy particle info from the GADGET file, with added coordinates \code{$galaxy_obs$z_obs},
#'   \code{$galaxy_obs$r_obs}, and \code{$galaxy_obs$vy_obs}. \code{$galaxy_obs$binn} gives the cell
#'   index in which each particle exists in the final cube. Depending on the \code{\link{sim_data}}
#'   file provided, the Age and Metallicity of particles will be included if SSP has been requested,
#'   Luminosity will be included per particle if SSP is not requested.}
#' \item{\code{$z}}{The projected redshift at which the observation is made.}
#' \item{\code{$sbin}}{The number of spatial bins.}
#' \item{\code{$sbinsize}}{The size of the spatial bins in kpc.}
#' \item{\code{$angular_size}}{The angular size of the galaxy in kpc/arcecond at the provided
#'  redshift.}
#' \item{\code{$vbin}}{The number of velocity bins.}
#' \item{\code{$vbinsize}}{The size of the velocity bins in km/s.}
#' \item{\code{$lsf_size}}{The gaussian standard deviation of the line spread function in km/s.}
#' \item{\code{$ap_region}}{The aperture region mask used to remove flux outside of the specified
#'  aperture.}
#'@examples
#' galaxy_data = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))
#' output = obs_data_prep(simdata      = galaxy_data,
#'                        z            = 0.1,
#'                        fov          = 15,
#'                        ap_shape     = "circular",
#'                        central_wvl  = 4800,
#'                        lsf_fwhm     = 2.65,
#'                        pixel_sscale = 0.5,
#'                        pixel_vscale = 1.04,
#'                        inc_deg      = 0)
#'

obs_data_prep = function(simdata, r200=200, z=0.05, fov=15, ap_shape="circular", central_wvl=4800,
                         lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04, inc_deg=70, align=FALSE){

  set.seed(42);
  sol_lum = 3.827e33;                                        # solar luminosity in erg s-1

  ang_size     = celestial::cosdistAngScale(z, ref="Planck") # angular size given z, kpc/"
  ap_size      = ang_size * fov                              # diameter size of the telescope, kpc
  sbin         = floor(fov / pixel_sscale)                   # number of spatial bins
  sbinsize     = ap_size / sbin                              # spatial bin size (kpc per bin)
  vbinsize     = (pixel_vscale / central_wvl) * (3e8 / 1e3)  # km/s per velocity bin
  lsf_size     = ((lsf_fwhm / central_wvl) * (3e8 / 1e3)) / (2 * sqrt(2*log(2))) # velocity uncertainty (sd)

  if (ap_shape == "circular"){
    ap_region = .circular_ap(sbin)
  }  # circular apperture mask
  if (ap_shape == "square"){
    ap_region = matrix(data = 1, ncol = sbin, nrow = sbin)
  }    # square apperture mask
  if (ap_shape == "hexagonal"){
    ap_region = .hexagonal_ap(sbin)
  } # hexagonal apperture mask

  result = grepl(paste(c("PartType2", "PartType3", "PartType4"), collapse = "|"), names(simdata))
  # finding the luminous matter within the simulation for imaging
  #SSP = FALSE # placeholder for SSP for loop, remains FALSE if no SSP in simdata

  if (any(result)) {
    present = which(result)
    if (length(present) == 1){
      galaxy_data = simdata[[present[1]]]$Part
      # if ("SSP" %in% names(simdata[[present[1]]])){
      #   SSP = which(grepl("SSP", names(simdata[[present[1]]])))}
    } else if (length(present) == 2){
      galaxy_data = rbind(simdata[[present[1]]]$Part, simdata[[present[2]]]$Part)
      # if (any("SSP" %in% c(names(simdata[[present[1]]]), names(simdata[[present[2]]])))){
      #   SSP = which(grepl("SSP", c(names(simdata[[present[1]]]), names(simdata[[present[2]]]))))}
    } else if (length(present) == 3){
      galaxy_data = rbind(simdata[[present[1]]]$Part, simdata[[present[2]]]$Part, simdata[[present[3]]]$Part)
      # if (any("SSP" %in% c(names(simdata[[present[1]]]), names(simdata[[present[2]]]), names(simdata[[present[3]]])))){
      #   SSP = which(grepl("SSP", c(names(simdata[[present[1]]]), names(simdata[[present[2]]]), names(simdata[[present[3]]]))))}
    }
  } else {
    cat("There are no particles representing luminous matter in this simulation (i.e. no stars, bulge or disc particles). \n")
    stop("LumPart Error") # if no stars, bulge or disc component, stop trying to build IFU cube
  }
  # concatenating the seperate ptype into a single particle data frame

  galaxy_data = cen_galaxy(galaxy_data) # centre data
  if (align){
    galaxy_data = .reorient_galaxy(galaxy_data) # reorient galaxies such that the semimajor axis lies along the x
  }
  galaxy_df   = obs_galaxy(galaxy_data, inc_deg * (pi / 180)) # extracting the position and LOS data
  galaxy_cdf  = galaxy_df[(galaxy_df$r < r200),] # removing particles beyond r200

  vbin = ceiling((max(abs(galaxy_cdf$vy_obs))*2) / vbinsize) # number of velocity bins
  if (vbin <= 2){vbin = 3}

  pixel_index  = seq(1,sbin*sbin*vbin, by=1)
  pixel_region  = array(ap_region, dim = c(sbin,sbin,vbin)) * pixel_index

  vseq = seq(-(vbin * vbinsize) / 2,
             (vbin * vbinsize) / 2, by=vbinsize) # velocity bin break positions

  sseq = seq(-(sbin * sbinsize) / 2,
             (sbin * sbinsize) / 2, by=sbinsize) # spatial bin break positions

  galaxy_cdf$binn = cut(galaxy_cdf$x, breaks=sseq, labels=F) +
    (sbin * cut(galaxy_cdf$z_obs, breaks=sseq, labels=F)) +
    (sbin^2 * cut(galaxy_cdf$vy_obs, breaks=vseq, labels=F)) - (sbin^2 + sbin) # assigning particles to positions in cube

  galaxy_cdf = galaxy_cdf[!is.na(galaxy_cdf$binn),]
  galaxy_cdf = galaxy_cdf[galaxy_cdf$binn %in% pixel_region,]

  if (length(present) == 1){
    if ("SSP" %in% names(simdata[[present[1]]])){
      galaxy_cdf$Metallicity = simdata[[present[1]]]$SSP[simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID,]$Metallicity
      galaxy_cdf$Age         = simdata[[present[1]]]$SSP[simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID,]$Age
    }
    else if ("Lum" %in% names(simdata[[present[1]]])){
      galaxy_cdf$Lum = simdata[[present[1]]]$Lum[as.integer(simdata[[present[1]]]$Part$ID) %in% as.integer(galaxy_cdf$ID)]
    }
  } else if (length(present) == 2){
      if (all("SSP" %in% c(names(simdata[[present[1]]]), names(simdata[[present[2]]])))){
        galaxy_cdf$Metallicity = c(simdata[[present[1]]]$SSP[simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID,]$Metallicity,
                                   simdata[[present[2]]]$SSP[simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID,]$Metallicity)
        galaxy_cdf$Age         = c(simdata[[present[1]]]$SSP[simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID,]$Age,
                                   simdata[[present[2]]]$SSP[simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID,]$Age)

      }
      else if (all("Lum" %in% c(names(simdata[[present[1]]]), names(simdata[[present[2]]])))){
        galaxy_cdf$Lum = c(simdata[[present[1]]]$Lum[as.integer(simdata[[present[1]]]$Part$ID) %in% as.integer(galaxy_cdf$ID)],
                           simdata[[present[2]]]$Lum[as.integer(simdata[[present[2]]]$Part$ID) %in% as.integer(galaxy_cdf$ID)])
      }
  } else if (length(present) == 3){
      if (all("SSP" %in% c(names(simdata[[present[1]]]), names(simdata[[present[2]]]), names(simdata[[present[3]]])))){
        galaxy_cdf$Metallicity = c(simdata[[present[1]]]$SSP[simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID,]$Metallicity,
                                   simdata[[present[2]]]$SSP[simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID,]$Metallicity,
                                   simdata[[present[3]]]$SSP[simdata[[present[3]]]$Part$ID %in% galaxy_cdf$ID,]$Metallicity)
        galaxy_cdf$Age         = c(simdata[[present[1]]]$SSP[simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID,]$Age,
                                   simdata[[present[2]]]$SSP[simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID,]$Age,
                                   simdata[[present[3]]]$SSP[simdata[[present[3]]]$Part$ID %in% galaxy_cdf$ID,]$Age)

      }
      else if (all("Lum" %in% c(names(simdata[[present[1]]]), names(simdata[[present[2]]]), names(simdata[[present[3]]])))){
        galaxy_cdf$Lum = c(simdata[[present[1]]]$Lum[as.integer(simdata[[present[1]]]$Part$ID) %in% as.integer(galaxy_cdf$ID)],
                           simdata[[present[2]]]$Lum[as.integer(simdata[[present[2]]]$Part$ID) %in% as.integer(galaxy_cdf$ID)],
                           simdata[[present[3]]]$Lum[as.integer(simdata[[present[3]]]$Part$ID) %in% as.integer(galaxy_cdf$ID)])
      }
  }
  # adding the relevant Lum/SSP info to the trimmed particles for later processing

  output = list("galaxy_obs"  = galaxy_cdf,
                "z"           = z,
                "sbin"        = sbin,
                "sbinsize"    = sbinsize,
                "angular_size"= ang_size,
                "vbin"        = vbin,
                "vbinsize"    = vbinsize,
                "lsf_size"    = lsf_size,
                "ap_region"   = ap_region)

  return(output)
}
