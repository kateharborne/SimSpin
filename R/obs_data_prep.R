# Kate Harborne (last edit - 29/11/2018)
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
#'@param filter When stellar particles and SEDs are provided, the filter used when converting
#'spectra into luminosity. Options include "r" for SDSS r filter and "g" for the SDSS g filter.
#'@return Returned is a list that contains:
#' \item{\code{$galaxy_obs}}{A data frame of the observed particle information (which contains the
#'   galaxy particle info from the GADGET file, with added coordinates \code{$galaxy_obs$z_obs},
#'   \code{$galaxy_obs$r_obs}, and \code{$galaxy_obs$vy_obs} and particle fluxes,
#'   \code{$galaxy_obs$flux}, in units of 1e-16 erg s-1 cm-2).}
#' \item{\code{$sbin}}{The number of spatial bins.}
#' \item{\code{$sbinsize}}{The size of the spatial bins in kpc.}
#' \item{\code{$angular_size}}{The angular size of the galaxy in kpc/arcecond at the provided
#'  redshift.}
#' \item{\code{$vbin}}{The number of velocity bins.}
#' \item{\code{$vbinsize}}{The size of the velocity bins in km/s.}
#' \item{\code{$lsf_size}}{The gaussian standard deviation of the line spread function in km/s.}
#' \item{\code{$appregion}}{The aperture region mask used to remove flux outside of the specified
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
                         lsf_fwhm=2.65, pixel_sscale=0.5, pixel_vscale=1.04, inc_deg=70, align=TRUE){

  set.seed(42);
  sol_lum = 3.827e33;                                        # solar luminosity in erg s-1

  ang_size     = celestial::cosdistAngScale(z, ref="Planck") # angular size given z, kpc/"
  lum_dist     = celestial::cosdistLumDist(z, ref="Planck")  # the luminosity distance, Mpc
  ap_size      = ang_size * fov                              # diameter size of the telescope, kpc
  sbin         = floor(fov / pixel_sscale)                   # bin sizes in the x- & y/z_obs- axes
  sbinsize     = ap_size / sbin                              # kpc per bin
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
  galaxy_df   = galaxy_df[(galaxy_df$r < r200),] # removing particles beyond r200

  if (ap_shape == "circular"){
    galaxy_cdf  = .circular_ap_cut(galaxy_df, ap_size)
  } # removing particles outside aperture
  if (ap_shape == "square"){
    galaxy_cdf = .square_ap_cut(galaxy_df, sbin, sbinsize)
  }
  if (ap_shape == "hexagonal"){
    galaxy_cdf = .hexagonal_ap_cut(galaxy_df, sbin, sbinsize)
  }

  vbin = ceiling((max(galaxy_cdf$vy_obs) - min(galaxy_cdf$vy_obs)) / vbinsize)

  sseq = seq(-(sbin * sbinsize) / 2,
                        (sbin * sbinsize) / 2, by=sbinsize)                   # spatial bin break positions
  # vseq = seq(-(vbin * vbinsize) / 2, (vbin * vbinsize) / 2, by=vbinsize)    # velocity bin break positions

  galaxy_cdf$binn = cut(galaxy_cdf$x, breaks=sseq, labels=F) +
    (sbin * cut(galaxy_cdf$z_obs, breaks=sseq, labels=F)) - sbin              # assigning particles to positions in aperture

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

  # If SSP is present, assign relevant Age/Metallicities. Else, assign luminoisty.
  #  if (any(as.logical(SSP))){
  #   if (length(SSP) == 1){
  #     inc_ID = which(simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID)
  #     galaxy_cdf$Age = simdata[[present[1]]]$SSP$Age[inc_ID] * 1e9
  #     galaxy_cdf$Metallicity = simdata[[present[1]]]$SSP$Metallicity[inc_ID]
  #   } else if (length(SSP) == 2){
  #     inc_ID_1 = which(simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID)
  #     inc_ID_2 = which(simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID)
  #     galaxy_cdf$Age = c(simdata[[present[1]]]$SSP$Age[inc_ID_1] * 1e9, simdata[[present[2]]]$SSP$Age[inc_ID_2] * 1e9)
  #     galaxy_cdf$Metallicity = c(simdata[[present[1]]]$SSP$Metallicity[inc_ID_1], simdata[[present[2]]]$SSP$Metallicity[inc_ID_2])
  #   } else if (length(SSP) == 3){
  #     inc_ID_1 = which(simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID)
  #     inc_ID_2 = which(simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID)
  #     inc_ID_3 = which(simdata[[present[3]]]$Part$ID %in% galaxy_cdf$ID)
  #     galaxy_cdf$Age = c(simdata[[present[1]]]$SSP$Age[inc_ID_1] * 1e9, simdata[[present[2]]]$SSP$Age[inc_ID_2] * 1e9, simdata[[present[3]]]$SSP$Age[inc_ID_3] * 1e9)
  #     galaxy_cdf$Metallicity = c(simdata[[present[1]]]$SSP$Metallicity[inc_ID_1], simdata[[present[2]]]$SSP$Metallicity[inc_ID_2], simdata[[present[3]]]$SSP$Metallicity[inc_ID_3])}
  # } else {
  #   if (length(present) == 1){
  #     inc_ID = which(simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID)
  #     galaxy_cdf$flux = (simdata[[present[1]]]$Lum[inc_ID] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2)
  #   } else if (length(present) == 2){
  #     inc_ID_1 = which(simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID)
  #     inc_ID_2 = which(simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID)
  #     galaxy_cdf$flux = c((simdata[[present[1]]]$Lum[inc_ID_1] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2),
  #                         (simdata[[present[2]]]$Lum[inc_ID_2] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2))
  #   } else if (length(present) == 3){
  #     inc_ID_1 = which(simdata[[present[1]]]$Part$ID %in% galaxy_cdf$ID)
  #     inc_ID_2 = which(simdata[[present[2]]]$Part$ID %in% galaxy_cdf$ID)
  #     inc_ID_3 = which(simdata[[present[3]]]$Part$ID %in% galaxy_cdf$ID)
  #     galaxy_cdf$flux = c((simdata[[present[1]]]$Lum[inc_ID_1] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2),
  #                         (simdata[[present[2]]]$Lum[inc_ID_2] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2),
  #                         (simdata[[present[3]]]$Lum[inc_ID_3] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2))
  #
  #   }
  # }

  # tic()
  # ID_grid = vector("list", sbin*sbin)
  # for (i in 1:length(galaxy_cdf$ID)){
  #   id = galaxy_cdf$binn[i]
  #   ID_grid[[id]] = c(galaxy_cdf$ID[i], ID_grid[[id]])
  # }
  # toc()

  # the number of velocity bins
  # galaxy_cdf$flux = rep(0, nrow(galaxy_cdf))
  #
  # if ("PartType2" %in% names(simdata)){                    # if there are disc particles in the data file,
  #   inc_discID = which(simdata$PartType2$Part$ID %in% galaxy_cdf$ID)
  #   disc_ids = simdata$PartType2$Part$ID[inc_discID]
  #   galaxy_cdf[(galaxy_cdf$ID %in% disc_ids),]$flux = (simdata$PartType2$Lum[inc_discID] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2)
  #   # calculating flux from disc particles fluxes in units of (1e-16 erg s-1 cm-2 arcsecond-2)
  # }
  #
  # if ("PartType3" %in% names(simdata)){                    # if there are bulge particles in the data file,
  #   inc_bulgeID = which(simdata$PartType3$Part$ID %in% galaxy_cdf$ID)
  #   bulge_ids = simdata$PartType3$Part$ID[inc_bulgeID]
  #   galaxy_cdf[(galaxy_cdf$ID %in% bulge_ids),]$flux = (simdata$PartType3$Lum[inc_bulgeID] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2)
  #   # calculating flux from bulge particles fluxes in units of (1e-16 erg s-1 cm-2 arcsecond-2)
  # }

  # if ("PartType4" %in% names(simdata)){                    # if there are stars in the data file,
  #
  #   inc_starID = which(simdata$PartType4$Part$ID %in% galaxy_cdf$ID)
  #   star_ids = simdata$PartType4$Part$ID[inc_starID]
  #
  #   if (length(names(simdata$PartType4)) == 3){              #  and if there are spectra associated,
  #     if (filter == "g"){tempfilt=list(ProSpect::filt_g_SDSS)} #  use ProSpect to get particle flux.
  #     if (filter == "r"){tempfilt=list(ProSpect::filt_r_SDSS)}
  #     star_flux = numeric(length = length(galaxy_cdf$ID))
  #     j = 1
  #     for (i in inc_starID){
  #       galaxy_cdf[(galaxy_cdf$ID %in% star_ids),]$flux[j] =
  #         ProSpect::photom_lum(wave = simdata$PartType4$Wav, lum = simdata$PartType4$Lum[i,], z = z,
  #                              outtype = "Jansky", filters = tempfilt, ref="Planck") # flux density in Janksy
  #       j = j+1
  #     }
  #   } else {
  #     galaxy_cdf[(galaxy_cdf$ID %in% star_ids),]$flux = (simdata$PartType4$Lum[inc_starID] * sol_lum * (pixel_sscale^2)) / (1e-16 * 4 * pi * (lum_dist * 3.086e24)^2)
  #   }
  # }

  output = list("galaxy_obs"  = galaxy_cdf,
                "sbin"        = sbin,
                "sbinsize"    = sbinsize,
                "angular_size"= ang_size,
                "vbin"        = vbin,
                "vbinsize"    = vbinsize,
                "lsf_size"    = lsf_size,
                "ap_region"   = ap_region,
                "pixel_sscale"= pixel_sscale,
                "d_L"         = lum_dist)

  return(output)
}
