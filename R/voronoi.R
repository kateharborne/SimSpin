# Kate Harborne 20/09/21
# Voronoi tessellation of pixels to meet some minimum particle number threshold


.make_bin = function(part_in_spaxel, each, bin_id){

  if (part_in_spaxel$N[each] >= observation$particle_limit){
    if(verbose){print(paste("Pixel", each, "meets limit - no binning."))}
    part_in_spaxel$bin_number[each] = bin_id
    bin_id = bin_id + 1
  } else {
    if(verbose){print(paste("Pixel", each, "does not meet limit - binning commencing."))}

    # comupte reff, where reff is the radius of a circle with the same area as the whole bin
    reff = sqrt(2 / pi) # need to join 2 pixels (min)

    unbinned_spaxels = part_in_spaxel[(each+1):total_pixels,]
    unbinned_spaxels$xbin_cen = unbinned_spaxels$xbin - part_in_spaxel$xbin[each]
    unbinned_spaxels$ybin_cen = unbinned_spaxels$ybin - part_in_spaxel$ybin[each]
    unbinned_spaxels$rad = sqrt(unbinned_spaxels$xbin_cen^2 + unbinned_spaxels$ybin_cen^2)

    possible_candidates = unbinned_spaxels[unbinned_spaxels$rad <= sqrt(2),] # meets topological criterion

    if(verbose){print(paste(length(possible_candidates$pixel_pos), "pixels meet topological criterion (connected pixels)."))}
    if (length(possible_candidates$pixel_pos) > 1){
      possible_candidates$roundness = (possible_candidates$rad / reff) - 1 # using the roundness parameter in Cappellari & Copin (2003)

      possible_candidates = possible_candidates[which(possible_candidates$roundness == min(possible_candidates$roundness) &
                                                        possible_candidates$roundness <= roundness_limit),]

      if(verbose){print(paste(length(possible_candidates$pixel_pos), "pixels meet morphological criterion (roundness)."))}

      if (length(possible_candidates$pixel_pos) > 1){
        possible_candidates$uniformity = abs((possible_candidates$N + part_in_spaxel$N[each]) - observation$particle_limit)
        join_bin = possible_candidates[which(possible_candidates$uniformity == min(possible_candidates$uniformity))[1],]
      }else{join_bin = possible_candidates}
    }else{join_bin = possible_candidates}

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

    # then remove the joined pixel from the array
    part_in_spaxel = part_in_spaxel[-id_join,]

    return(part_in_spaxel)

  }

}

voronoi = function(part_in_spaxel, observation, roundness_limit = 0.3){

  total_pixels = length(part_in_spaxel$pixel_pos)

  # Following the algorithm set out in Cappellari & Copin (2003):
  # (i) Begin at the highest valued pixel in the image.

  # order by number of particles (max first)
  data.table::setorder(part_in_spaxel, -N)

  # x-y position of pixel
  part_in_spaxel$xbin = part_in_spaxel$pixel_pos %% observation$sbin
  part_in_spaxel$ybin = part_in_spaxel$pixel_pos %/% observation$sbin

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

      # comupte reff, where reff is the radius of a circle with the same area as the whole bin
      reff = sqrt(2 / pi) # need to join 2 pixels (min)

      unbinned_spaxels = part_in_spaxel[(each+1):total_pixels,]
      unbinned_spaxels$xbin_cen = unbinned_spaxels$xbin - part_in_spaxel$xbin[each]
      unbinned_spaxels$ybin_cen = unbinned_spaxels$ybin - part_in_spaxel$ybin[each]
      unbinned_spaxels$rad = sqrt(unbinned_spaxels$xbin_cen^2 + unbinned_spaxels$ybin_cen^2)

      possible_candidates = unbinned_spaxels[unbinned_spaxels$rad <= sqrt(2),] # meets topological criterion

      if(verbose){print(paste(length(possible_candidates$pixel_pos), "pixels meet topological criterion (connected pixels)."))}
      if (length(possible_candidates$pixel_pos) > 1){
        possible_candidates$roundness = (possible_candidates$rad / reff) - 1 # using the roundness parameter in Cappellari & Copin (2003)

        possible_candidates = possible_candidates[which(possible_candidates$roundness == min(possible_candidates$roundness) &
                                                          possible_candidates$roundness <= roundness_limit),]

        if(verbose){print(paste(length(possible_candidates$pixel_pos), "pixels meet morphological criterion (roundness)."))}

        if (length(possible_candidates$pixel_pos) > 1){
          possible_candidates$uniformity = abs((possible_candidates$N + part_in_spaxel$N[each]) - observation$particle_limit)
          join_bin = possible_candidates[which(possible_candidates$uniformity == min(possible_candidates$uniformity))[1],]
        }else{join_bin = possible_candidates}
      }else{join_bin = possible_candidates}

      # Having identified a bin to join, let's add the particles to that bin
      part_in_spaxel$bin_number[each] = bin_id
      # which particles are in the full bin?
      part_in_bin = c(part_in_spaxel$val[[each]], part_in_spaxel$val[[which(part_in_spaxel$pixel_pos == join_bin$pixel_pos)]])
      part_in_spaxel$val[[each]] = part_in_bin


      part_in_spaxel$bin_number[each] = bin_id
      part_in_spaxel$bin_number[part_in_spaxel$pixel_pos == join_bin$pixel_pos] = bin_id

      if (sum(part_in_spaxel$N[part_in_spaxel$bin_number == bin_id]) >= observation$particle_limit){
        bin_id = bin_id + 1

      } else {
        stop()
      }

    }

  }

  if (sum(part_in_spaxel$N[part_in_spaxel$bin_number == bin_id]) >= observation$particle_limit){
    bin_id = bin_id + 1

  } else {
    if(verbose){print(paste("More pixels need adding to bin to meet vorbin_limit!"))}
  }

}
