# Author: Kate Harborne
# Date: 26/10/2020
# Title: Function for generating a data cube from the observation

build_datacube = function(simspin_file, telescope, observing_strategy){

  observation = observation(telescope = telescope, observing_strategy = observing_strategy)

  galaxy_data = data.table::transpose(fst::read_fst(simspin_file, from = 1, to = 7)[,-1]) #read in position and velocity data
  colnames(galaxy_data) = c("ID", "x", "y", "z", "vx", "vy", "vz")

  galaxy_data = obs_galaxy(part_data = galaxy_data, inc_rad = observation$inc_rad) # projecting the galaxy to given inclination

  galaxy_data$pixel_pos = cut(galaxy_data$x, breaks=observation$sbin_seq, labels=F) +
    (observation$sbin * cut(galaxy_data$z_obs, breaks=observation$sbin_seq, labels=F)) - (observation$sbin) # assigning particles to positions in cube

  galaxy_data = galaxy_data[!is.na(galaxy_data$pixel_pos),] # removing any particles that fall outside the sbin aperture

  galaxy_data = galaxy_data[galaxy_data$pixel_pos %in% observation$aperture_region,] # trimming particles that lie outside the aperture of the telescope

  for (i in range(galaxy_data$pixel_pos)){
    particle_IDs = galaxy_data$ID[galaxy_data$pixel_pos == i]
    spectra = fst::read_fst(simspin_file, columns = particle_IDs, from = 8)
  }

}
