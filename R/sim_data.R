# Kate Harborne (last edit - 23/04/2018)
#'Prepare simulation data for SimSpin
#'

sim_data = function(filename, ptype=NA, SSP=NA, m2l_disc=1, m2l_bulge=1, m2l_star = 1){

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
    disc_part = data.frame("ID"        = as.integer(paste0("2", disc_ID)),
                           "x"         = hdf5r::readDataSet(galaxy_file[["PartType2/x"]]),
                           "y"         = hdf5r::readDataSet(galaxy_file[["PartType2/y"]]),
                           "z"         = hdf5r::readDataSet(galaxy_file[["PartType2/z"]]),
                           "vx"        = hdf5r::readDataSet(galaxy_file[["PartType2/vx"]]),
                           "vy"        = hdf5r::readDataSet(galaxy_file[["PartType2/vy"]]),
                           "vz"        = hdf5r::readDataSet(galaxy_file[["PartType2/vz"]]),
                           "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType2/Mass"]]))
    disc_lum = (disc_part$Mass * 1e10) / m2l_disc
    PartType2 = list("Part" = disc_part, "Lum" = disc_lum)

  }

  if ("PartType3" %in% ptype){ # if bulge particles are requested
    bulge_n = length(hdf5r::readDataSet(galaxy_file[["PartType3/x"]]))
    bulge_ID = formatC(1:bulge_n, width = floor(log10(bulge_n)) + 1, format = "d", flag = "0")
    bulge_part = data.frame("ID"        = as.integer(paste0("3", bulge_ID)),
                            "x"         = hdf5r::readDataSet(galaxy_file[["PartType3/x"]]),
                            "y"         = hdf5r::readDataSet(galaxy_file[["PartType3/y"]]),
                            "z"         = hdf5r::readDataSet(galaxy_file[["PartType3/z"]]),
                            "vx"        = hdf5r::readDataSet(galaxy_file[["PartType3/vx"]]),
                            "vy"        = hdf5r::readDataSet(galaxy_file[["PartType3/vy"]]),
                            "vz"        = hdf5r::readDataSet(galaxy_file[["PartType3/vz"]]),
                            "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType3/Mass"]]))
    bulge_lum = (bulge_part$Mass * 1e10) / m2l_bulge
    PartType3 = list("Part" = bulge_part, "Lum" = bulge_lum)

  }

  if ("PartType4" %in% ptype){ # if gas particles are requested
    star_n = length(hdf5r::readDataSet(galaxy_file[["PartType4/x"]]))
    star_ID = formatC(1:star_n, width = floor(log10(star_n)) + 1, format = "d", flag = "0")
    star_part = data.frame("ID"        = as.integer(paste0("4", star_ID)),
                           "x"         = hdf5r::readDataSet(galaxy_file[["PartType4/x"]]),
                           "y"         = hdf5r::readDataSet(galaxy_file[["PartType4/y"]]),
                           "z"         = hdf5r::readDataSet(galaxy_file[["PartType4/z"]]),
                           "vx"        = hdf5r::readDataSet(galaxy_file[["PartType4/vx"]]),
                           "vy"        = hdf5r::readDataSet(galaxy_file[["PartType4/vy"]]),
                           "vz"        = hdf5r::readDataSet(galaxy_file[["PartType4/vz"]]),
                           "Mass"      = hdf5r::readDataSet(galaxy_file[["PartType4/Mass"]]))

    if (!is.na(SSP)){
      f = hdf5r::h5file(SSP, mode = "r")
      star_lum = hdf5r::readDataSet(f[["PartType4/Luminosity"]])
      star_wave = hdf5r::readDataSet(f[["PartType4/Wavelength"]])
      hdf5r::h5close(f)
      PartType4 = list("Part" = star_part, "Lum" = star_lum, "Wav" = star_wave)
    } else {
      star_lum = ((star_part$Mass * 1e10) / m2l_star)
      PartType4 = list("Part" = star_part, "Lum" = star_lum)
    }
  }

  hdf5r::h5close(galaxy_file)                                 # close the snapshot data file

  galaxy_data = setNames(lapply(ls(pattern="PartType*"), function(x) get(x)), ls(pattern="PartType*"))

  return(galaxy_data)
}
