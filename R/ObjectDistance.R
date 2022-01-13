




setClass("ObjectDistance",
         slots = c(
           distance = "numeric",
           unit = "character"
           ),
         prototype = list(
           distance = NA_real_,
           unit = NA_character_
         )
        )

setValidity("ObjectDistance", function(object){

  if (length(object@distance) != length(object@unit)) {
    "@distance and @unit must be the same length"
    # check that passed parameters are the same length
  } else if (object@distance < 0) {
    "@distance cannot be negative"
    # check that the @distance provided is never 0
  } else if (stringr::str_to_upper(object@unit) != "Z" & stringr::str_to_upper(object@unit) != "MPC"
             & stringr::str_to_upper(object@unit) != "KPC/ARCSEC") {
    "@units must be specified as 'z', 'Mpc' or 'kpc/arcsec'"
    # check that a valid unit definition is provided
  } else if (stringr::str_to_upper(object@unit) == "Z" & object@distance == 0) {
    "@distance cannot be zero when @unit is redshift, z"
    # if distance is given as a redshift, @distance cannot be 0
  } else if (stringr::str_to_upper(object@unit) == "KPC/ARCSEC" & object@distance == 0) {
    "@distance cannot be zero when @unit is angular, kpc/arcsec"
    # if distance is specified as an angular size, @distance cannot be 0
  } else {
    TRUE
  }

})

# Setting generics and methods for modifying/printing the distance of ObjectDistance
setGeneric("distance", function(x) standardGeneric("distance"))
setMethod("distance", "ObjectDistance", function(x) x@distance)

# Assessors for units of ObjectDistance
setGeneric("distance<-", function(x, value) standardGeneric("distance<-"))
setMethod("distance<-", "ObjectDistance", function(x, value) {
  x@distance <- value
  validObject(x)
  x
})

# Setting generics and methods for modifying/printing the units of ObjectDistance
setGeneric("unit", function(x) standardGeneric("unit"))
setMethod("unit", "ObjectDistance", function(x) x@unit)

# Assessors for units of ObjectDistance
setGeneric("unit<-", function(x, value) standardGeneric("unit<-"))
setMethod("unit<-", "ObjectDistance", function(x, value) {
  x@unit <- value
  validObject(x)
  x
})



ObjectDistance <- function(distance, unit){

  distance = as.double(distance)

  new("ObjectDistance", distance = distance, unit = unit)

}
