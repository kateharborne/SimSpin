# Kate Harborne 20/09/21
# Voronoi tessellation of pixels to meet some minimum particle number threshold


.accretion = function(x, y, count, target_count, pixel_size){

  n = length(x)
  classe = logical(length=n)
  good = numeric(length=n)


}

voronoi = function(part_in_spaxel, observation){

  # Following the algorithm set out in Cappellari & Copin (2003):
  # (i) Begin at the highest valued pixel in the image.

  # number of particles per pixel
  part_in_spaxel$count = sapply(X = part_in_spaxel$val, FUN = length)
  # order by number of particles (max first)
  data.table::setorder(part_in_spaxel, -count)
  # x-y position of pixel
  part_in_spaxel$xbin = part_in_spaxel$spaxel_ID %% observation$sbin
  part_in_spaxel$ybin = part_in_spaxel$spaxel_ID %/% observation$sbin

  # setting "bin_number" value
  part_in_spaxel$bin_number = 0
  bin_id = 1
  total_pixels = length(part_in_spaxel$spaxel_ID)
  for (each in 1:total_pixels){

    if (part_in_spaxel$count[each] >= observation$particle_limit){
      if(verbose){print(paste("Pixel", each, "meets limit - no binning."))}
      part_in_spaxel$bin_number[each] = bin_id
      bin_id = bin_id + 1
    } else {
      if(verbose){print(paste("Pixel", each, "does not meet limit - binning commencing."))}

      reff = (2 / pi) # need to join 2 pixels (min)
      unbinned_spaxels = part_in_spaxel[(each+1):total_pixels,]
      unbinned_spaxels$xbin = unbinned_spaxels$xbin - part_in_spaxel$xbin[each]
      unbinned_spaxels$ybin = unbinned_spaxels$ybin - part_in_spaxel$ybin[each]
      unbinned_spaxels$rad2 = unbinned_spaxels$xbin^2 + unbinned_spaxels$ybin^2

      possible_candidates = unbinned_spaxels[unbinned_spaxels$rad2 <= 2,] # meets topological criterion
      if(verbose){print(paste(length(possible_candidates$spaxel_ID), "pixels meet topological criterion (connected pixels)."))}
      if (length(possible_candidates$spaxel_ID) > 1){
        possible_candidates$roundness = (possible_candidates$rad2 / reff) - 1
        possible_candidates = possible_candidates[which(possible_candidates$roundness == min(possible_candidates$roundness)),]
        if(verbose){print(paste(length(possible_candidates$spaxel_ID), "pixels meet morphological criterion (roundness)."))}
        if (length(possible_candidates$spaxel_ID) > 1){
          possible_candidates$uniformity = abs((possible_candidates$count + part_in_spaxel$count[each]) - observation$particle_limit)
          join_bin = possible_candidates[which(possible_candidates$uniformity == min(possible_candidates$uniformity))[1],]
        }else{join_bin = possible_candidates}
      }else{join_bin = possible_candidates}

      part_in_spaxel$bin_number[each] = bin_id
      part_in_spaxel$bin_number[part_in_spaxel$spaxel_ID == join_bin$spaxel_ID] = bin_id

      if (sum(part_in_spaxel$count[part_in_spaxel$bin_number == bin_id]) > observation$particle_limit){
        bin_id = bin_id + 1
      } else {
        stop()
      }

    }

  }

}
