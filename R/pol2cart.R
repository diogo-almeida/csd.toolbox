#' @title Translates polar to cartesian coordinates
#' 
#' @description Adaptation of the MATLAB function of the same name.
#'
#' @param Theta Theta coordinates in polar space (in radians).
#'
#' @param Phi Phi coordinates in polar space.
#'
#' @return List containing the corresponding x- and y-locations in cartesian 
#'   coordinates
#'
#' @export
#'
#' @seealso \code{\link{sph2cart}}
#'
#' @examples
#' theta <- 129.254
#' phi   <- 29.833
#' pol2cart(theta, phi)
#'
#' thetas <- c(129.254, 90, 50.746)
#' phis   <- c(29.833, 45, 29.833)
#' pol2cart(thetas, phis)
pol2cart <- function(Theta, Phi){
   x1 <- Phi * cos(Theta)
   y1 <- Phi * sin(Theta)
   XY <- list(x = x1, y = y1)
   return(XY)
}
