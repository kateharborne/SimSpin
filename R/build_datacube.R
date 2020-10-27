# Author: Kate Harborne
# Date: 26/10/2020
# Title: Function for generating a data cube from the observation

build_datacube = function(simspin_file, telescope, observing_strategy, verbose = F, write_fits = F, output_location, ...){

  if (verbose){cat("Computing observation parameters... \n")}
  observation = observation(telescope = telescope, observing_strategy = observing_strategy)

  galaxy_data = data.table::transpose(fst::read_fst(simspin_file, from = 1, to = 7)[,-1]) #read in position and velocity data
  colnames(galaxy_data) = c("ID", "x", "y", "z", "vx", "vy", "vz")

  galaxy_data = obs_galaxy(part_data = galaxy_data, inc_rad = observation$inc_rad) # projecting the galaxy to given inclination

  if (verbose){cat("Assigning particles to spaxels... \n")}
  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z_obs, breaks=observation$sbin_seq, labels=F)) - (observation$sbin) # assigning particles to positions in cube

  galaxy_data = galaxy_data[!is.na(galaxy_data$pixel_pos),] # removing any particles that fall outside the sbin aperture

  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$aperture_region,] # trimming particles that lie outside the aperture of the telescope

  wavelength  = fst::read_fst(simspin_file, columns = "V1", from = 8)[,1]

  spectra = matrix(data = NA, ncol = observation$wave_bin, nrow = observation$sbin^2)

  if (verbose){cat("Generating spectra per spaxel... \n")}
  for (i in sort(unique(galaxy_data$pixel_pos))){ # computing the spectra at each occupied spatial pixel position
     particle_IDs = paste0("V", (galaxy_data$ID[galaxy_data$pixel_pos == i] + 1)) # IDs are indexed one out from the column numbers
     intrinsic_spectra = fst::read_fst(simspin_file, columns = particle_IDs, from = 8)
     velocity_los = galaxy_data$vy_obs[galaxy_data$pixel_pos == i]
     wave = matrix(data = rep(wavelength, length(particle_IDs)), nrow = length(particle_IDs), byrow=T)
     wave_shift = ((velocity_los / observation$c) * wave) + wave # using doppler formula to compute the shift in wavelengths cause by velocity
     spectra[i,] = .interpolate_spectra(shifted_wave = wave_shift, spectra = intrinsic_spectra, wave_seq = observation$wave_seq)
     if (verbose){cat(i, "... ", sep = "")}
  }

  cube = array(data = spectra, dim = c(observation$sbin, observation$sbin, observation$wave_bin))
  if (verbose){cat("Done! \n")}

  if (write_fits){
    if (verbose){cat("Writing FITS... \n")}
    if (missing(output_location)){
      out_file_name = stringr::str_remove(simspin_file, ".fst")
      output_location = paste(out_file_name, "_inc", observation$inc_deg, "deg_seeing", observation$psf_fwhm,"fwhm.FITS", sep="")
    }

    crpixn = dim(cube)%/%2
    crvaln = c(observation$sbin_seq[crpixn[1]]/observation$ang_size/3600,
               observation$sbin_seq[crpixn[2]]/observation$ang_size/3600,
               observation$wave_seq[crpixn[3]])
    cdeltn = c(observation$sbin_size/observation$ang_size/3600,
               observation$sbin_size/observation$ang_size/3600,
               observation$wave_res)
    ctypen = c("RA---TAN", "DEC--TAN", "AWAV")
    cunitn = c("deg", "deg", "Angstrom")
    axDat  = data.frame("crpix" = crpixn, "crval" = crvaln, "cdelt" = cdeltn,
                        "len" = dim(cube), "ctype" = ctypen, "cunit" = cunitn)

    FITSio::writeFITSim(X = cube, file = output_location, crpixn = crpixn,
                        crvaln = crvaln, cdeltn = cdeltn, ctypen = ctypen,
                        cunitn = cunitn, axDat = axDat)

    message("FITS file written to: ", output_location)
  }

  output = list("spectral_cube" = cube,
                "observation"   = observation)

  return(output)
}

.interpolate_spectra = function(shifted_wave, spectra, wave_seq){

  shifted_spectra = vector(mode = "list", length = dim(shifted_wave)[1])
  for(j in 1:dim(shifted_wave)[1]){
    shifted_spectra[[j]] = approx(x = shifted_wave[j,], y = spectra[,j], xout = wave_seq, rule=1)[2]
  }

  output = matrix(unlist(shifted_spectra, use.names=FALSE), nrow = dim(shifted_wave)[1], byrow = TRUE)
  output = colSums(output)
  return(output)
}


