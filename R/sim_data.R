# Kate Harborne (last edit - 03/12/2018)
#'Prepare simulation data for SimSpin
#'
#'The purpose of this function is to take a SimSpin simulation file (an HDF5 file with particle
#' information subset in particle types, "PartTypeX") and output a data.frame that is easy for
#' processing throughout the rest of the package. This reduces the number of times that you need to
#' load the simulation into the package, which is the most time-consuming process. This also allows
#' you to combine the particle information with spectral information output from other codes.
#'
#'@param filename The SimSpin HDF5 file containing the particle information of the galaxy.
#'@param ptype The particle type/types to be extracted - NA (default) gives all particles in the
#' simulation, 0 - gas, 1 - dark matter, 2 - disc, 3 - bulge, 4 - stars, 5 - boundary.
#'@param SSP Boolean specifying whether you wish the flux of the particle be calculated based
#' on the supplied Age/Metallicity information or scaled by the mass. Default is FALSE. If set
#' to TRUE and no Age/Metallicity information is present in the file, function will error.
#'@param m2l_disc The mass-to-light ratio of the disc component in solar units.
#'@param m2l_bulge The mass-to-light ratio of the bulge component in solar units.
#'@param m2l_star If no SSP file is specified and stellar particles exist in the file, the
#' mass-to-light ratio of the stellar component in solar units.
#'
#'@return A list of data.frames containing the particle information (\code{$Part}) for each
#' particle type requested from the simulation. Each data frame contains the position (x, y, z)
#' and velocity (vx, vy, vz) information, along with the ID and mass of each particle. Also
#' associated with each PartType element in the list is an array of luminosities, \code{$Lum}. If
#' SSP information has been supplied, instead of Luminosities, a data frame \code{$SSP} will be
#' included in the list which contains the Metallicity and Ages of each particle such that spectra
#' can be constructed at later part in the process.
#'@examples
#' output = sim_data(system.file("extdata", 'SimSpin_example.hdf5', package="SimSpin"))

sim_data = function(filename, ptype=NA, SSP=FALSE, m2l_disc=1, m2l_bulge=1, m2l_star=1){

  galaxy_file = hdf5r::h5file(filename, mode = "r")           # reading in the snapshot data
  ppart = substring(hdf5r::list.groups(galaxy_file), 1)
  if (is.na(ptype[1])){
    ptype = ppart # if NA is chosen, extract all present particle types
  } else {
    ptype = paste("PartType", as.character(ptype), sep="")
    # else convert ptype into a comparible format
    if (!all(ptype %in% ppart)){ # if any requested particle types are not within the file
      cat("Particles of", paste(ptype[!ptype %in% ppart], collapse = ","), "are missing in this model. \n")
      stop("Npart Error")
    }
  }

  if ("PartType0" %in% ptype){ # if gas particles are requested
    gas_n = length(hdf5r::readDataSet(galaxy_file[["PartType0/x"]]))
    gas_ID = formatC(1:gas_n, width = floor(log10(gas_n)) + 1, format = "d", flag = "0")
    gas_part = data.frame("ID"        = as.integer(paste0("9", gas_ID)),
                          "x"         = hdf5r::readDataSet(galaxy_file[["PartType0/x"]]),
                          "y"         = hdf5r::readDataSet(galaxy_file[["PartType0/y"]]),
                          "z"         = hdf5r::readDataSet(galaxy_file[["PartType0/z"]]),
                          "vx"        = hdf5r::readDataSet(galaxy_file[["PartType0/vx"]]),
                          "vy"        = hdf5r::readDataSet(galaxy_file[["PartType0/vy"]]),
                          "vz"        = hdf5r::readDataSet(galaxy_file[["PartType0/vz"]]),
                          "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType0/Mass"]]))
    PartType0 = list("Part" = gas_part)
  }

  if ("PartType1" %in% ptype){ # if DM particles are requested
    DM_n = length(hdf5r::readDataSet(galaxy_file[["PartType1/x"]]))
    DM_ID = formatC(1:DM_n, width = floor(log10(DM_n)) + 1, format = "d", flag = "0")
    DM_part = data.frame("ID"        = as.integer(paste0("1", DM_ID)),
                         "x"         = hdf5r::readDataSet(galaxy_file[["PartType1/x"]]),
                         "y"         = hdf5r::readDataSet(galaxy_file[["PartType1/y"]]),
                         "z"         = hdf5r::readDataSet(galaxy_file[["PartType1/z"]]),
                         "vx"        = hdf5r::readDataSet(galaxy_file[["PartType1/vx"]]),
                         "vy"        = hdf5r::readDataSet(galaxy_file[["PartType1/vy"]]),
                         "vz"        = hdf5r::readDataSet(galaxy_file[["PartType1/vz"]]),
                         "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType1/Mass"]]))
    PartType1 = list("Part" = DM_part)
  }

  if ("PartType2" %in% ptype){ # if disc particles are requested
    disc_n = length(hdf5r::readDataSet(galaxy_file[["PartType2/x"]]))
    disc_ID = formatC(1:disc_n, width = floor(log10(disc_n)) + 1, format = "d", flag = "0")

    if (length(hdf5r::list.datasets(galaxy_file[["PartType2"]])) == 7 && SSP){
      cat("SSP requested, but no Age/Metallicity information contained within supplied file,", paste(filename), ". \n",
          "Please set SSP=FALSE, or provide additional particle information. \n")
      stop("SSP Error")        # if SSP is requested without the required info, error
    } else {
      disc_part = data.frame("ID"        = as.integer(paste0("2", disc_ID)),
                             "x"         = hdf5r::readDataSet(galaxy_file[["PartType2/x"]]),
                             "y"         = hdf5r::readDataSet(galaxy_file[["PartType2/y"]]),
                             "z"         = hdf5r::readDataSet(galaxy_file[["PartType2/z"]]),
                             "vx"        = hdf5r::readDataSet(galaxy_file[["PartType2/vx"]]),
                             "vy"        = hdf5r::readDataSet(galaxy_file[["PartType2/vy"]]),
                             "vz"        = hdf5r::readDataSet(galaxy_file[["PartType2/vz"]]),
                             "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType2/Mass"]]))
    }

     if (length(hdf5r::list.datasets(galaxy_file[["PartType2"]])) == 7 && isFALSE(SSP)){
      disc_lum = (disc_part$Mass * 1e10) / m2l_disc
      PartType2 = list("Part" = disc_part, "Lum" = disc_lum)

    } else if (length(hdf5r::list.datasets(galaxy_file[["PartType2"]])) == 9 && isFALSE(SSP)){
      disc_lum = (disc_part$Mass * 1e10) / m2l_disc
      PartType2 = list("Part" = disc_part, "Lum" = disc_lum)

    } else if (length(hdf5r::list.datasets(galaxy_file[["PartType2"]])) == 9 && SSP){

      if ("Age" %in% substring(hdf5r::list.datasets(galaxy_file[["PartType2"]]), 1)){
        disc_SSP = data.frame("Metallicity" = hdf5r::readDataSet(galaxy_file[["PartType2/Metallicity"]]),
                              "Age"         = hdf5r::readDataSet(galaxy_file[["PartType2/Age"]]))

        # if 9 datasets are provided and SSP is true, check if "Age" is already supplied
      } else { # if not, calculate the "Age" using celestial (performed in parallel)
        numCores = parallel::detectCores()
        sft = hdf5r::readDataSet(galaxy_file[["PartType2/StellarFormationTime"]])
        age = as.numeric(parallel::mclapply(sft, .SFTtoAge, mc.cores = numCores))

        disc_SSP = data.frame("Metallicity" = hdf5r::readDataSet(galaxy_file[["PartType2/Metallicity"]]),
                              "Age" = age)
      }

      PartType2 = list("Part" = disc_part, "SSP" = disc_SSP)

    }
  }

  if ("PartType3" %in% ptype){ # if bulge particles are requested
    bulge_n = length(hdf5r::readDataSet(galaxy_file[["PartType3/x"]]))
    bulge_ID = formatC(1:bulge_n, width = floor(log10(bulge_n)) + 1, format = "d", flag = "0")
    if (length(hdf5r::list.datasets(galaxy_file[["PartType3"]])) == 7 && SSP){
      cat("SSP requested, but no Age/Metallicity information contained within supplied file,", paste(filename), ". \n",
          "Please set SSP=FALSE, or provide additional particle information. \n")
      stop("SSP Error")        # if SSP is requested without the required info, error
    } else {
      bulge_part = data.frame("ID"        = as.integer(paste0("3", bulge_ID)),
                              "x"         = hdf5r::readDataSet(galaxy_file[["PartType3/x"]]),
                              "y"         = hdf5r::readDataSet(galaxy_file[["PartType3/y"]]),
                              "z"         = hdf5r::readDataSet(galaxy_file[["PartType3/z"]]),
                              "vx"        = hdf5r::readDataSet(galaxy_file[["PartType3/vx"]]),
                              "vy"        = hdf5r::readDataSet(galaxy_file[["PartType3/vy"]]),
                              "vz"        = hdf5r::readDataSet(galaxy_file[["PartType3/vz"]]),
                              "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType3/Mass"]]))
    }

    if (length(hdf5r::list.datasets(galaxy_file[["PartType3"]])) == 7 && isFALSE(SSP)){
      bulge_lum = (bulge_part$Mass * 1e10) / m2l_bulge
      PartType3 = list("Part" = bulge_part, "Lum" = bulge_lum)

    } else if (length(hdf5r::list.datasets(galaxy_file[["PartType3"]])) == 9 && isFALSE(SSP)){
      bulge_lum = (bulge_part$Mass * 1e10) / m2l_bulge
      PartType3 = list("Part" = bulge_part, "Lum" = bulge_lum)

    } else if (length(hdf5r::list.datasets(galaxy_file[["PartType3"]])) == 9 && SSP){

      if ("Age" %in% substring(hdf5r::list.datasets(galaxy_file[["PartType3"]]), 1)){
        bulge_SSP = data.frame("Metallicity" = hdf5r::readDataSet(galaxy_file[["PartType3/Metallicity"]]),
                               "Age"         = hdf5r::readDataSet(galaxy_file[["PartType3/Age"]]))

        # if 9 datasets are provided and SSP is true, check if "Age" is already supplied
      } else { # if not, calculate the "Age" using celestial (performed in parallel)
        numCores = parallel::detectCores()
        sft = hdf5r::readDataSet(galaxy_file[["PartType3/StellarFormationTime"]])
        age = as.numeric(parallel::mclapply(sft, .SFTtoAge, mc.cores = numCores))

        bulge_SSP = data.frame("Metallicity" = hdf5r::readDataSet(galaxy_file[["PartType3/Metallicity"]]),
                               "Age" = age)
      }

      PartType3 = list("Part" = bulge_part, "SSP" = bulge_SSP)

    }
  }

  if ("PartType4" %in% ptype){ # if stellar particles are requested

    star_n = length(hdf5r::readDataSet(galaxy_file[["PartType4/x"]]))
    star_ID = formatC(1:star_n, width = floor(log10(star_n)) + 1, format = "d", flag = "0")

    if (length(hdf5r::list.datasets(galaxy_file[["PartType4"]])) == 7 && SSP){
      cat("SSP requested, but no Age/Metallicity information contained within supplied file,", paste(filename), ". \n",
      "Please set SSP=FALSE, or provide additional particle information. \n")
      stop("SSP Error")        # if SSP is requested without the required info, error
    } else {
      star_part = data.frame("ID"        = as.integer(paste0("4", star_ID)),
                             "x"         = hdf5r::readDataSet(galaxy_file[["PartType4/x"]]),
                             "y"         = hdf5r::readDataSet(galaxy_file[["PartType4/y"]]),
                             "z"         = hdf5r::readDataSet(galaxy_file[["PartType4/z"]]),
                             "vx"        = hdf5r::readDataSet(galaxy_file[["PartType4/vx"]]),
                             "vy"        = hdf5r::readDataSet(galaxy_file[["PartType4/vy"]]),
                             "vz"        = hdf5r::readDataSet(galaxy_file[["PartType4/vz"]]),
                             "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType4/Mass"]]))

    }

    if (length(hdf5r::list.datasets(galaxy_file[["PartType4"]])) == 7 && isFALSE(SSP)){

      star_lum = ((star_part$Mass * 1e10) / m2l_star)
      # if no Age/Metallicity is provided, use M2L ratio to calculate luminosity
      PartType4 = list("Part" = star_part, "Lum" = star_lum)

    } else if (length(hdf5r::list.datasets(galaxy_file[["PartType4"]])) == 10 && isFALSE(SSP)){

      star_lum = ((star_part$Mass * 1e10) / m2l_star)
      # if Age/Metallicity is provided, but SSP is false, use M2L ratio to calculate luminosity
      PartType4 = list("Part" = star_part, "Lum" = star_lum)

    } else if (length(hdf5r::list.datasets(galaxy_file[["PartType4"]])) == 10 && SSP){

      if ("Age" %in% substring(hdf5r::list.datasets(galaxy_file[["PartType4"]]), 1)){

        star_SSP = data.frame("Metallicity" = hdf5r::readDataSet(galaxy_file[["PartType4/Metallicity"]]),
                              "Age"         = hdf5r::readDataSet(galaxy_file[["PartType4/Age"]]),
                              "InitialMass" = hdf5r::readDataSet(galaxy_file[["PartType4/InitialMass"]]))

        # if 9 datasets are provided and SSP is true, check if "Age" is already supplied
      } else { # if not, calculate the "Age" using celestial (performed in parallel)

        numCores = parallel::detectCores()
        sft = hdf5r::readDataSet(galaxy_file[["PartType4/StellarFormationTime"]])
        age = as.numeric(parallel::mclapply(sft, .SFTtoAge, mc.cores = numCores))

        star_SSP = data.frame("Metallicity" = hdf5r::readDataSet(galaxy_file[["PartType4/Metallicity"]]),
                              "Age" = age,
                              "InitialMass" = hdf5r::readDataSet(galaxy_file[["PartType4/InitialMass"]]))

      }

      PartType4 = list("Part" = star_part, "SSP" = star_SSP)

    }
  }

  hdf5r::h5close(galaxy_file)                                 # close the snapshot data file

  galaxy_data = setNames(lapply(ls(pattern="PartType*"), function(x) get(x)), ls(pattern="PartType*"))

  return(galaxy_data)
}
