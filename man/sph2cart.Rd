\name{sph2cart}
\alias{sph2cart}
\title{Translates spherical coordinates to cartesian coordinates}
\usage{
  sph2cart(Theta, Phi, R)
}
\arguments{
  \item{Theta}{Azimuth in spherical coordinates.}

  \item{Phi}{Elevation in spherical coordinates.}

  \item{R}{R in spherical coordinates.}
}
\value{
  List containing the corresponding x-, y- and z-locations
  in cartesian coordinates
}
\description{
  Adaptation of the MATLAB function of the same name.
}
\examples{
theta <- 129.254
phi   <- 29.833
r     <- 1
sph2cart(theta, phi, r)

thetas <- c(129.254, 90, 50.746)
phis   <- c(29.833, 45, 29.833)
r      <- rep(1,3)
sph2cart(thetas, phis, r)
}
\seealso{
  \code{\link{pol2cart}}
}

