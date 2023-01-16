# Kate Harborne 20/09/21
# Voronoi tessellation of pixels to meet some minimum particle number threshold

.make_bin = function(part_in_spaxel, each, bin_id, verbose){

    # define some flags for whether or not a bin can accrete any more pixels
    topo_criteria = F
    roundness_criteria = F
    uniform_criteria = F

    # comupte reff, where reff is the radius of a circle with the same area as the whole bin
    reff = sqrt((part_in_spaxel$area[each] + 1) / pi) # need to join bin to at least one pixel (min)

    unbinned_spaxels = part_in_spaxel[(each+1):total_pixels,]
    unbinned_spaxels$xbin_cen = unbinned_spaxels$xbin - part_in_spaxel$xbin[each]
    unbinned_spaxels$ybin_cen = unbinned_spaxels$ybin - part_in_spaxel$ybin[each]
    unbinned_spaxels$rad = sqrt(unbinned_spaxels$xbin_cen^2 + unbinned_spaxels$ybin_cen^2)

    possible_candidates = unbinned_spaxels[unbinned_spaxels$rad <= sqrt(2),] # do we meet the topological criterion?

    if (length(possible_candidates$pixel_pos) >= 1){
      topo_criteria = T
      if(verbose){print(paste(length(possible_candidates$pixel_pos), "pixels meet topological criterion (connected pixels)."))}
      } else {topo_criteria = F}

    if (topo_criteria){ # if we meet the topological criterion, test for roundness

      possible_candidates$roundness = (possible_candidates$rad / reff) - 1 # using the roundness parameter in Cappellari & Copin (2003)

      possible_candidates = possible_candidates[which(possible_candidates$roundness == min(possible_candidates$roundness) &
                                                        possible_candidates$roundness <= roundness_limit),]

      if (length(possible_candidates$pixel_pos) >= 1){ # do we meet roundness?
        roundness_criteria = T
        if(verbose){print(paste(length(possible_candidates$pixel_pos), "pixels meet morphological criterion (roundness)."))}
      } else {roundness_criteria = F}

      if (roundness_criteria){ # if we meet the roundness criteria, test for uniformity

        possible_candidates$uniformity = abs(1 - ((possible_candidates$N + part_in_spaxel$N[each])/observation$particle_limit))

        possible_candidates = possible_candidates[which(possible_candidates$uniformity <= 0.2),]

        if (length(possible_candidates$pixel_pos) >= 1){ # do we meet uniformity?
          uniform_criteria = T
          join_bin = possible_candidates[which(possible_candidates$uniformity == min(possible_candidates$uniformity))[1],]
          # if we also meet the uniformity criteria, add the bin that gets us closest to the target
        } else {uniform_criteria = F}

      }

    }

    if (topo_criteria && roundness_criteria && uniform_criteria){ # if we find a suitable pixel to add to the bin

    id_join = which(part_in_spaxel$pixel_pos == join_bin$pixel_pos) # joining pixel

    # Having identified a pixel to join, let's add the particles to that bin
    part_in_spaxel$bin_number[each] = bin_id

    # which particles are in the bin?
    part_in_bin = c(part_in_spaxel$val[[each]], part_in_spaxel$val[[id_join]])
    part_in_spaxel$val[[each]] = part_in_bin # particles in bin updated

    xbin_cen = median(c(part_in_spaxel$xbin[each], part_in_spaxel$xbin[id_join]))
    ybin_cen = median(c(part_in_spaxel$ybin[each], part_in_spaxel$ybin[id_join]))

    # updating to the centre of the new bin
    part_in_spaxel$xbin[each] = xbin_cen
    part_in_spaxel$ybin[each] = ybin_cen

    # updating the area of the bin
    new_area = part_in_spaxel$area[each] + part_in_spaxel$area[id_join]
    part_in_spaxel$area[each] = new_area

    # then remove the joined pixel from the array
    part_in_spaxel = part_in_spaxel[-id_join,]

    } else { # flag this pixel as unable to join to another
      part_in_spaxel$bin_number[each] = NA
    }

    return(part_in_spaxel)

}

voronoi = function(part_in_spaxel, observation, roundness_limit = 0.3, uniform_limit = 0.8){

  total_pixels = length(part_in_spaxel$pixel_pos)

  # Following the algorithm set out in Cappellari & Copin (2003):
  # (i) Begin at the highest valued pixel in the image.

  # order by number of particles (max first)
  data.table::setorder(part_in_spaxel, -N)

  # x-y position of pixel
  part_in_spaxel$xbin = part_in_spaxel$pixel_pos %% observation$sbin
  part_in_spaxel$ybin = part_in_spaxel$pixel_pos %/% observation$sbin
  part_in_spaxel$area = 1

  # setting "bin_number" value
  part_in_spaxel$bin_number = 0
  bin_id = 1

  for (each in 1:total_pixels){

    if (part_in_spaxel$N[each] >= observation$particle_limit){
      if(verbose){print(paste("Pixel", each, "meets limit - no binning."))}
      part_in_spaxel$bin_number[each] = bin_id
      bin_id = bin_id + 1
    } else {
      if(verbose){print(paste("Pixel", each, "does not meet limit - binning commencing."))}

      part_in_spaxel = .make_bin(part_in_spaxel, each, bin_id, verbose)

      while (part_in_spaxel$N[each] < (observation$particle_limit*uniform_limit) &&
             !is.na(part_in_spaxel$bin_number[each])){
        if(verbose){print(paste("More pixels need adding to bin to meet vorbin_limit!"))}
        part_in_spaxel = .make_bin(part_in_spaxel, each, bin_id, verbose)
      }

      if(verbose && !is.na(part_in_spaxel$bin_number[each])){print(paste("DONE adding pixels to bin to meet vorbin_limit!"))}
      if(verbose && is.na(part_in_spaxel$bin_number[each])){print(paste("No pixels available. Flagging pixel and moving on."))}
      bin_id = bin_id + 1

    }

    readline(prompt=paste("Each is ", each, ". Press [enter] to continue"))

  }


}
