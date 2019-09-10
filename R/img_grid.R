

img_grid = function(obs_data, z, filter="g", threshold){

  galIDs = as.data.frame(obs_data$galaxy_obs$binn) # spaxel that each particle sits in
  sbin = obs_data$sbin # dimensions of the final image = sbin*sbin
  image_grid = lapply(seq(1,sbin*sbin), function(x) which(galIDs == x)) # list of particles that exist in each spaxel
  lengths = lapply(image_grid, length) # number of particles in each spaxel
  flux = numeric(sbin*sbin) # flux in each spaxel
  if (filter=="g"){tempfilt=list(ProSpect::filt_g_SDSS)} # loading the filter
  if (filter=="r"){tempfilt=list(ProSpect::filt_r_SDSS)}

  for (cell in 1:(sbin*sbin)){ # in each pixel
    if (lengths[[cell]] > 1){ # if there are particles in the cell
      temp = data.frame("Age" = obs_data$galaxy_obs$Age[image_grid[[cell]]] * 1e9,
                        "Mass" = obs_data$galaxy_obs$Mass[image_grid[[cell]]] * 1e10,
                        "cumMass" = rep(NA, length(obs_data$galaxy_obs$Mass[image_grid[[cell]]])),
                        "Metallicity" = obs_data$galaxy_obs$Metallicity[image_grid[[cell]]])
      temp = temp[order(temp$Age),] # fill a data frame and order by age

      if (length(unique(temp$Age)) != length(temp$Age)){
        # if there are any particles with the same age
        nu = as.integer(which(table(temp$Age) > 1))
        for (n in 1:length(nu)){
          vals = which(temp$Age == temp$Age[nu[n]])
          temp$Mass[vals[1]] = sum(temp$Mass[vals])
          temp$Metallicity[vals[1]] = sum(temp$Metallicity[vals])
          vals = vals[-1]
          temp = temp[-vals,]
        } # sum together all masses and metallicities such that we have a single particle at that age
      }

      temp$cumMass = cumsum(temp$Mass) # cumulative sum of ages
      tmp_SFHfunc = approxfun(x = c(temp$Age), y = c(0,diff(temp$cumMass)/diff(temp$Age)), yleft=0, yright=0)
      # describe star formation rate in M_sol/yr as a function of age
      tmp_Zfunc = approxfun(x = c(temp$Age), y = c(temp$Metallicity), yleft=0, yright=0)
      # describe metallicity as a function of age
      tmp_SFH = ProSpect::SFHfunc(massfunc = tmp_SFHfunc, Z = tmp_Zfunc, z=z, filters = tempfilt, ref="Planck")
      # generate a spectrum
      flux[cell] = ProSpect::photom_lum(wave = tmp_SFH$wave_lum, lum = tmp_SFH$lum_unatten, z=z, outtype="Jansky", filters=tempfilt, ref="Planck")
      # convert to flux
    }
  }

  threshold_flux  = ProSpect::magAB2Jansky(threshold)
  flux[flux < threshold_flux] = 0 # remove any cells that aren't above the threshold
  flux_1 = matrix(data = flux, nrow=sbin, ncol=sbin, byrow=TRUE) # NEED TO CHECK WHICH TO USE!!!
  flux_2 = matrix(data = flux, nrow=sbin, ncol=sbin, byrow=FALSE)
  return(flux)

}
