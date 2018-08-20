# Kate Harborne 09/08/18
#'Kinematic analysis of the LOSVD for simulation data.
#'
#'The purpose of this function is to fit the parameters of the line-of-sight velocity distribution
#' (LOSVD) and return the gauss-hermite coefficients, h$_3$ and h$_4$.
#'
#'@param filename The Gadget file containing the particle information of the galaxy to be analysed.
#'@param inc_deg The inclination at which to observe the galaxy in degrees.
#'@param rmax The maximum radius to consider when constructing the LOSVD.
#'@param ptype The particle type/types to be extracted - \code{NA} (default) gives all particles in
#' the simulation, 1 - gas, 2 - dark matter, 3 - disc, 4 - bulge, 5 - stars, 6 - boundary.
#'@param los_upper The upper boundary for the L-BFGS-B fit.
#'@param los_lower The lower boundary for the L-BFGS-B fit.
#'@param plot_out Logical specifying whether or not you would like the function to output a plot. 
#' Default is \code{TRUE}.
#'@return A list that contains the particle data of the simulation (\code{$part_data} and 
#' \code{$head_data} for a particle data.frame and header file respectively), a data frame 
#' containing the projected data, and the best fit h$_3$ and h$_4$ coefficients.
#' 

sim_losvd = function(filename, inc_deg, rmax=200, ptype=NA, los_upper=c(1,1), los_lower=c(-1,-1), plot_out = TRUE){
  
  galaxy_data = snapshot::snapread(filename)
  
  galaxy_data$part$part_type = rep(0, nrow(galaxy_data$part))
                                                           # add a "particle type" column
  p = seq(1,6)                                             # all possible particle values
  ppart = cumsum(galaxy_data$head$Npart[which(galaxy_data$head$Nall[p] != 0)])
                                                           # number of present particles
  for (i in 1:length(ppart)){
    if (i == 1){
      galaxy_data$part[1:as.integer(ppart[i]),]$part_type =  which(galaxy_data$head$Nall[p] != 0)[i]
    } else {
      galaxy_data$part[as.integer(ppart[i-1]+1):as.integer(ppart[i]),]$part_type = which(galaxy_data$head$Nall[p] != 0)[i]
    }
  }
                                                           # labelling the data with particle types
  
  galaxy_df = obs_galaxy(part_data = galaxy_data$part, centre = TRUE, inc_rad = inc_deg * (pi / 180))
  
  galaxy_df = galaxy_df[galaxy_df$r < rmax,]               # removing particles beyond rmax
  
  losvd_lower = los_lower                                 # lower limit for fit
  losvd_upper = los_upper                                 # upper limit for fit
  
  galaxy_losvd = optim(.gh_sim, par = c(0, 0),
                       v_los = galaxy_df$vy_obs, m_v = mean(galaxy_df$vy_obs), s_v = sd(galaxy_df$vy_obs),
                       method='L-BFGS-B', lower=losvd_lower, upper=losvd_upper, hessian = TRUE)
                                                           # losvd fitting function
  
  d = approxfun(density(galaxy_df$vy_obs))
  dy = d(galaxy_df$vy_obs)
  l_v = .losvd(galaxy_df$vy_obs, m_v = mean(galaxy_df$vy_obs), s_v = sd(galaxy_df$vy_obs), h3 = galaxy_losvd$par[1], h4 = galaxy_losvd$par[2])
  if (isTRUE(plot_out)){
    h = magicaxis::maghist(galaxy_df$vy_obs, breaks = 50, plot = FALSE, verbose = FALSE)
    magicaxis::maghist(galaxy_df$vy_obs, breaks = 50, scale = 1/max(h$counts), freq=TRUE,
                       ylab ="", xlab=expression("v"[LOS]*", km s"^-1), main = "", #ylim=c(0, max(l_v)),
                       verbose = FALSE)
    points(galaxy_df$vy_obs, l_v/max(l_v), col="red", pch=".", cex = 3)
    legend("topright", inset=c(0,0.01), c(parse(text=sprintf('h[3] == %f',galaxy_losvd$par[1])),parse(text=sprintf('h[4] == %f',galaxy_losvd$par[2]))) , text.col = "black", bty="n")
    
  }
  
  return(list("part_data" = galaxy_data$part, 
              "head_data" = galaxy_data$head, 
              "projected_profile" = galaxy_df, 
              "h3" = galaxy_losvd$par[1], 
              "h4" = galaxy_losvd$par[2]))
  
}
  
