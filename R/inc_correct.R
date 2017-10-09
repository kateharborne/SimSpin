# Kate Harborne (last edit - 14/09/2017)
#'Inclination correction for \eqn{\lambda}_R.
#'
#'The purpose of this function is to take the \eqn{\lambda}_R value observed via the \code{obs_lambda()} function and correct
#'it for inclination, as done for true IFU observations. You supply both the true and observed inclinations, as calculated
#'by the \code{obs_inclination()} function. The true anisotropy can also be supplied, as calculated from the simulation data
#'in \code{sim_analysis()}. Else, the inclination correction will be calculated for 4 anisotropy values between 0 and 0.6.
#'
#'@param inc_deg The true inclination at which the galaxy was observed (as specified in \code{obs_data_prep()}).
#'@param galaxy_class The Hubble type of the galaxy with options \code{c("S0","Sa", "Sb", "Sc", "Sd")}.
#'@param axis_ratio A data frame containing the semi-major and semi-minor axes lengths for the observed galaxy, as given by ifu_img().
#'@param obs_lambdar The observed \eqn{\lambda}_R value as produced by the \code{obs_lambda()} function.
#'@param v_aniso \emph{Optional} parameter that specifies the velocity anisotropy at which the inclination correction is made.
#' If missing, the correction will be applied for \code{v_aniso = c(0, 0.2, 0.4, 0.6)}.
#'@return Returns a list containing the observed inclination (\code{$obs_inc}) and a matrix containing the \eqn{\lambda}_R values (\code{$obs_lambda})
#' for a true inclination correction (\code{"true_inc_correct"}) and the observed inclination correction (\code{"obs_inc_correct"}) for each velocity
#' anisotropy.
#'@examples
#' \dontrun{
#' images      = ifu_img()
#' lambda_data = obs_lambda()
#'
#' inc_correct(inc_deg      = 0,
#'             galaxy_class = "S0",
#'             axis_ratio   = images$axis_ratio,
#'             obs_lambdar  = lambda_data$obs_lambdar)
#'
#' inc_correct(inc_deg      = 0,
#'             galaxy_class = "S0",
#'             axis_ratio   = images$axis_ratio,
#'             obs_lambdar  = lambda_data$obs_lambdar,
#'             v_aniso      = 0.1)
#'
#' inc_correct(inc_deg      = 0,
#'             galaxy_class = "S0",
#'             axis_ratio   = images$axis_ratio,
#'             obs_lambdar  = lambda_data$obs_lambdar,
#'             v_aniso     = c(0, 0.1, 0.2, 0.3))
#' }
#'

inc_correct = function(inc_deg, galaxy_class, axis_ratio, obs_lambdar, v_aniso = c(0, 0.2, 0.4, 0.6)){
  inc_rad     = inc_deg * (pi / 180)
  if (galaxy_class == "early"){
    fac = 0.6
  } else if (galaxy_class == "late"){
      fac = 0.2
  } else {
      cat("galaxy_class is not recognised. Please select from \"early\" or \"late\" and run again.")
      stop("galaxy_class Error")
    }
  obs_inc_rad = acos(sqrt(((axis_ratio$b / axis_ratio$a)^2 - fac^2) / (1 - fac^2))) # inclination from axis ratio as in 2016_cortese
  c_true      = sin(inc_rad) / (sqrt(1 - v_aniso * (cos(inc_rad)^2)))
  c_obs       = sin(obs_inc_rad) / (sqrt(1 - v_aniso * (cos(obs_inc_rad)^2)))
  lr_true     = (obs_lambdar / c_true) * (1 / (sqrt(1 - (obs_lambdar^2) * (1 - 1 / (c_true^2)))))
  lr_obs      = (obs_lambdar / c_obs) * (1 / (sqrt(1 - (obs_lambdar^2) * (1 - 1 / (c_obs^2)))))
  out1        = obs_inc_rad * (180 / pi)
  out2        = matrix(data = rbind(lr_true, lr_obs), nrow = 2, ncol = length(v_aniso), dimnames = list(c("true_inc_correct", "obs_inc_correct"), v_aniso))
  output      = list("obs_inc" = out1, "obs_lambda" = out2)
  return(output)
}
