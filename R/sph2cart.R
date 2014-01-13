#' @title Translates spherical coordinates to cartesian coordinates
#' 
#' @description Adaptation of the MATLAB function of the same name.
#'
#' @param Theta Azimuth in spherical coordinates.
#'
#' @param Phi Elevation in spherical coordinates.
#'
#' @param R R in spherical coordinates.
#'
#' @return List containing the corresponding x-, y- and z-locations in cartesian 
#'   coordinates
#'
#' @export
#'
#' @seealso \code{\link{pol2cart}}
#'
#' @examples
#' theta <- 129.254
#' phi   <- 29.833
#' r     <- 1
#' sph2cart(theta, phi, r)
#'
#' thetas <- c(129.254, 90, 50.746)
#' phis   <- c(29.833, 45, 29.833)
#' r      <- rep(1,3)
#' sph2cart(thetas, phis, r)
sph2cart <- function(Theta, Phi, R){
   x1  <- R * cos(Phi) * cos(Theta)
   y1  <- R * cos(Phi) * sin(Theta)
   z1  <- R * sin(Phi)
   XYZ <- list(x = x1, y = y1, z = z1)
   return(XYZ)
}
