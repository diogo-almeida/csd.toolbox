\name{pol2cart}
\alias{pol2cart}
\title{Translates polar to cartesian coordinates}
\usage{
  pol2cart(Theta, Phi)
}
\arguments{
  \item{Theta}{Theta coordinates in polar space (in
  radians).}

  \item{Phi}{Phi coordinates in polar space.}
}
\value{
  List containing the corresponding x- and y-locations in
  cartesian coordinates
}
\description{
  Adaptation of the MATLAB function of the same name.
}
\examples{
theta <- 129.254
phi   <- 29.833
pol2cart(theta, phi)

thetas <- c(129.254, 90, 50.746)
phis   <- c(29.833, 45, 29.833)
pol2cart(thetas, phis)
}
\seealso{
  \code{\link{sph2cart}}
}

