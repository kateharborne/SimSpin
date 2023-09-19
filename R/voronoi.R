# Kate Harborne 20/09/21
# Voronoi tessellation of pixels to meet some minimum particle number threshold

.make_bin = function(vorbin_spaxel, bin_id, verbose=F, index = NULL){

    # define some flags for whether or not a bin can accrete any more pixels
    topo_criteria = F
    roundness_criteria = F
    uniform_criteria = F

    if (is.null(index)){ # if this is the first call of make bin per round, start with all unbinned
      data.table::setorder(vorbin_spaxel, -N)
      index_to_bin = which(vorbin_spaxel$pixel_pos == vorbin_spaxel$pixel_pos[vorbin_spaxel$bin_number == 0 &
                                                                                !is.na(vorbin_spaxel$bin_number)][1])
      unbinned_spaxels = vorbin_spaxel[vorbin_spaxel$bin_number == 0,]
      unbinned_spaxels = unbinned_spaxels[-which(unbinned_spaxels$pixel_pos == vorbin_spaxel$pixel_pos[index_to_bin]),]

    } else { # if this is a repeat, consider the old index as the bin to grow
      index_to_bin = index
      unbinned_spaxels = vorbin_spaxel[vorbin_spaxel$bin_number == 0,]
    }

    # comupte reff, where reff is the radius of a circle with the same area as the whole bin
    reff = sqrt((vorbin_spaxel$area[index_to_bin] + 1) / pi) # need to join bin to at least one pixel (min)

    unbinned_spaxels$xbin_cen = unbinned_spaxels$xbin - vorbin_spaxel$xbin[index_to_bin]
    unbinned_spaxels$ybin_cen = unbinned_spaxels$ybin - vorbin_spaxel$ybin[index_to_bin]
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

        possible_candidates$uniformity = abs(1 - ((possible_candidates$N + vorbin_spaxel$N[index_to_bin])/observation$particle_limit))

        #possible_candidates = possible_candidates[which(possible_candidates$uniformity <= 0.2),]

        if (length(possible_candidates$pixel_pos) >= 1){ # do we meet uniformity?
          uniform_criteria = T
          join_bin = possible_candidates[which(possible_candidates$uniformity == min(possible_candidates$uniformity))[1],]
          # if we also meet the uniformity criteria, add the bin that gets us closest to the target
        } else {uniform_criteria = F}

      }

    }

    if (topo_criteria && roundness_criteria && uniform_criteria){ # if we find a suitable pixel to add to the bin

    id_join = which(vorbin_spaxel$pixel_pos == join_bin$pixel_pos) # joining pixel

    # Having identified a pixel to join, let's add the particles to that bin
    vorbin_spaxel$bin_number[index_to_bin] = bin_id

    # which particles are in the bin?
    part_in_bin = c(vorbin_spaxel$val[[index_to_bin]], vorbin_spaxel$val[[id_join]])
    vorbin_spaxel$val[[index_to_bin]] = part_in_bin # particles in bin updated

    # record ids of merged pixels
    vorbin_spaxel$joined_ids[index_to_bin][[1]] =
      c(vorbin_spaxel$joined_ids[index_to_bin][[1]], vorbin_spaxel$joined_ids[id_join][[1]])

    # update the new number of particles per bin
    vorbin_spaxel$N[index_to_bin] = sum(vorbin_spaxel$N[index_to_bin], vorbin_spaxel$N[id_join])

    xbin_cen = median(c(vorbin_spaxel$xbin[index_to_bin], vorbin_spaxel$xbin[id_join]))
    ybin_cen = median(c(vorbin_spaxel$ybin[index_to_bin], vorbin_spaxel$ybin[id_join]))

    # updating to the centre of the new bin
    vorbin_spaxel$xbin[index_to_bin] = xbin_cen
    vorbin_spaxel$ybin[index_to_bin] = ybin_cen

    # updating the area of the bin
    new_area = vorbin_spaxel$area[index_to_bin] + vorbin_spaxel$area[id_join]
    vorbin_spaxel$area[index_to_bin] = new_area

    # then remove the joined pixel in the array
    vorbin_spaxel = vorbin_spaxel[-id_join,]

    } else { # flag this pixel as unable to join to another
      vorbin_spaxel$bin_number[index_to_bin] = NA
    }

    return(list("vorbin"=vorbin_spaxel, "index" = index_to_bin))

}

voronoi = function(part_in_spaxel, observation, roundness_limit = 0.3, uniform_limit = 0.8){

  # Following the algorithm set out in Cappellari & Copin (2003):
  # (i) Begin at the highest valued pixel in the image.

  # order by number of particles (max first)
  data.table::setorder(part_in_spaxel, -N)

  # setting "bin_number" value
  part_in_spaxel$bin_number = 0

  vorbin_spaxels = part_in_spaxel
  # x-y position of pixel
  vorbin_spaxels$xbin = vorbin_spaxels$pixel_pos %% observation$sbin
  vorbin_spaxels$ybin = vorbin_spaxels$pixel_pos %/% observation$sbin
  vorbin_spaxels$area = 1
  vorbin_spaxels$joined_ids = list()

  for (pixel in 1:nrow(vorbin_spaxels)){
    vorbin_spaxels$joined_ids[pixel][[1]] = c(vorbin_spaxels$pixel_pos[pixel])
  }

  bin_id = 1
  each = 1

  # while pixels have not yet been binned, run algorithm
  while (nrow(vorbin_spaxels[vorbin_spaxels$bin_number == 0]) > 0){

    if (part_in_spaxel$N[each] >= observation$particle_limit){

      if(verbose){print(paste("Pixel", each, "meets limit - no binning."))}
      part_in_spaxel$bin_number[each] = bin_id
      vorbin_spaxels$bin_number[each] = bin_id

      bin_id = bin_id + 1
      each = each + 1

    } else {
      if(verbose){print(paste("Pixel", each, "does not meet limit - binning commencing."))}

      binned_spaxels = .make_bin(vorbin_spaxels, bin_id, verbose, index=NULL)
      vorbin_spaxels = binned_spaxels$vorbin
      index = binned_spaxels$index

      while (vorbin_spaxels$N[index] < (observation$particle_limit*uniform_limit) &&
             !is.na(vorbin_spaxels$bin_number[index])){
        if(verbose){print(paste("More pixels need adding to bin to meet vorbin_limit!"))}
        binned_spaxels = .make_bin(vorbin_spaxels, bin_id, verbose, index = index)
        vorbin_spaxels = binned_spaxels$vorbin
        index = binned_spaxels$index
      }

      if(!is.na(vorbin_spaxels$bin_number[index])){
        if(verbose){print(paste("DONE adding pixels to bin to meet vorbin_limit!"))}
        part_in_spaxel$bin_number[part_in_spaxel$pixel_pos %in% vorbin_spaxels$joined_ids[index][[1]]] =
          vorbin_spaxels$bin_number[index]
      }

      if(is.na(vorbin_spaxels$bin_number[index])){
        if(verbose){print(paste("No pixels available. Flagging pixel and moving on."))}
        # part_in_spaxel$bin_number[part_in_spaxel$pixel_pos %in% vorbin_spaxels$joined_ids[[index]]] =
        #   vorbin_spaxels$bin_number[index]
      }

      bin_id = bin_id + 1
      each = each + 1
    }

  }

  # iterate through pixels unable to be binned and add them to the nearest bin
  while (any(is.na(vorbin_spaxels$bin_number))){

    if(verbose){
      print(paste0(length(vorbin_spaxels$bin_number[is.na(vorbin_spaxels$bin_number)]), " unbinned pixels..."))
    }

    # take a single pixel that has not yet been binned
    unbinned_spaxel = vorbin_spaxels[is.na(vorbin_spaxels$bin_number),][1]
    bin_index = which(vorbin_spaxels$pixel_pos == unbinned_spaxel$pixel_pos)

    #compute distance between these and all binned pixels
    binned_spaxels = vorbin_spaxels[!is.na(vorbin_spaxels$bin_number) & vorbin_spaxels$bin_number > 0,]
    binned_spaxels$xbin_cen = binned_spaxels$xbin - unbinned_spaxel$xbin
    binned_spaxels$ybin_cen = binned_spaxels$ybin - unbinned_spaxel$ybin
    binned_spaxels$rad = sqrt(binned_spaxels$xbin_cen^2 + binned_spaxels$ybin_cen^2)

    join_to_bin = binned_spaxels$bin_number[which(binned_spaxels$rad == min(binned_spaxels$rad))]
    join_index = which(vorbin_spaxels$bin_number %in% join_to_bin)

    if (length(join_to_bin) > 1){ # if more than 1 close bin, minimise uniformity
      uniformity = abs(1 - ((binned_spaxels$N[join_index] + unbinned_spaxel$N)/observation$particle_limit))
      n = which(uniformity == min(uniformity))[1]
      join_index = join_index[n]
    } else {n = 1}

    # add unbinned to the closest bin
    vorbin_spaxels$joined_ids[join_index][[1]] =
      c(vorbin_spaxels$joined_ids[join_index][[1]], unbinned_spaxel$joined_ids[[1]])

    vorbin_spaxels$N[join_index] = sum(vorbin_spaxels$N[join_index],
                                       unbinned_spaxel$N)

    vorbin_spaxels$xbin[join_index] = median(vorbin_spaxels$xbin[join_index], unbinned_spaxel$xbin)
    vorbin_spaxels$ybin[join_index] = median(vorbin_spaxels$ybin[join_index], unbinned_spaxel$ybin)

    vorbin_spaxels$area[join_index] = sum(vorbin_spaxels$area[join_index], unbinned_spaxel$area)

    # update bin number of joined bin
    vorbin_spaxels = vorbin_spaxels[-bin_index,]

    # update the bin in part_in_spaxel
    part_in_spaxel$bin_number[part_in_spaxel$pixel_pos %in% vorbin_spaxels$joined_ids[join_index][[1]]] = join_to_bin[n]

  }

  # re-jigging the output part_in_spaxel
  unique_bins = unique(part_in_spaxel$bin_number)
  bn = length(unique_bins)
  new_part_in_spaxel = data.table::data.table("pixel_pos" = vector(mode="list", length=bn),
                                              "val" = vector(mode="list", length=bn),
                                              "N" = numeric(length=bn))

  for (bin in 1:bn){
    new_part_in_spaxel$pixel_pos[bin][[1]] = c(part_in_spaxel$pixel_pos[part_in_spaxel$bin_number == unique_bins[bin]])
    new_part_in_spaxel$val[bin][[1]] = c(unlist(part_in_spaxel$val[part_in_spaxel$bin_number == unique_bins[bin]]))
    new_part_in_spaxel$N[bin] = sum(part_in_spaxel$N[part_in_spaxel$bin_number == unique_bins[bin]])
  }

  return(new_part_in_spaxel)

}
