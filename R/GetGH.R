#' @title Computes the transformation matrices G and H for CSD transformation.
#' 
#' @description Computes the transformation matrices G and H to be used in the 
#'   CSD transformation.
#'
#' @param M Montage list (output of ExtractMontage function), containing the
#'   spherical and cartesian coordinates of the electrode montage to be used in
#'   the CSD transformation.
#'
#' @param m Value of the m-constant to be used in the spherical spline 
#'   interpolation. Lower values allow more flexible splines than higher values.
#'   The range of possible values goes from 2 to 10, 4 being the default. 
#'   Changing this value changes the properties of the spatial filter, so you 
#'   need to know what you are doing when changing this value.
#'
#' @return List containing the two transformation matrices G and H. Matrix G is
#'   used in computing the spherical spline interpolation of the surface 
#'   potentials, whereas matrix H is used in the spherical spline interpolation 
#'   that leads to the current source density estimate of the data.
#'
#' @note 
#'   This function requires computing Legendre polynomials. The original
#'   MATLAB code used the function \code{legendre}. The R base lacks a similar 
#'   function, but package \code{\link{orthopolynom}} implements an equivalent 
#'   solution, even though the user interface is substantially different
#'   from the MATLAB original function.
#' 
#'   IMPORTANT: The Legendre polynomial results from MATLAB and R are not 
#'   exactly the same all the time. The reasons at this point are unclear to me.
#'   However, despite not giving exactly the same results, tests based on the 
#'   tutorial data of the original MATLAB CSD Toolbox indicate that the MATLAB 
#'   output and the R output are equal up to the 6th and the 10th decimal places
#'   for G and H respectively, so the differences seem to be negligible. 
#'
#' @export
#'
#' @seealso \code{\link{ExtractMontage}}; \code{\link{CalculateCSD}}; 
#'   \code{\link{orthopolynom}}
#'
#' @examples
#'
#' montage <- ExtractMontage(colnames(NR.C66.trr))
#  gh      <- GetGH(montage)
#  str(gh)
GetGH <- function(M, m = 4) {
  require(orthopolynom)
  
  possible.m.values <- 2:10
  if (!(m %in% possible.m.values)) {
    stop(
      sprintf("Invalid m = %d (use 2, 3, 4, ..., 10)", m)
    )
  }

  ThetaRad <- (2 * pi * M$theta) / 360  # convert Theta and Phi to radians ...
  PhiRad   <- (2 * pi * M$phi) / 360    # ... and Cartesian coordinates ...
  XYZ <- sph2cart(ThetaRad, PhiRad, 1)  # ... for optimal resolution
  X   <- XYZ$x
  Y   <- XYZ$y
  Z   <- XYZ$z

  nElec <- length(M$electrode.labels)   # determine size of EEG montage
  # initialize interelectrode matrix...
  EF <- array(numeric(nElec * nElec), dim = c(nElec, nElec))
  # ...and compute all cosine distances
  for (i in 1:nElec) {
    for (j in 1:nElec) {
      EF[i, j] <- 1 - (((X[i] - X[j])^2 + (Y[i] - Y[j])^2 + (Z[i] - Z[j])^2 ) / 2)
    }
  }
  cat(sprintf("Spline flexibility:  m = %d\n", m))
  N <- 50  # set N iterations
  G <- array(numeric(nElec * nElec), dim = c(nElec,nElec)) # claim memory for 
  H <- array(numeric(nElec * nElec), dim = c(nElec,nElec)) # G- and H-matrices
  # intialize progress bar
  cat(sprintf("%d iterations for %d sites [", N, nElec)) 
  for (i in 1:nElec) {
    for (j in 1:nElec) {
      # compute Legendre polynomial
      # R Note (2009-11-14): this requires the R package "orthopolynom"
      #                      this part of the code changed quite a bit from the
      #                      original Matlab code.
      p.coefs <- legendre.polynomials(N, normalized = F)
      p <- as.numeric(polynomial.values(p.coefs, EF[i, j]))
      P <- p[2:(N + 1)]
      g = 0          # Initialize g and h in order to
      h = 0          # compute h- and g-functions
      if (j == 1) {
        cat('*')     # show progress
      }          
      for (n in 1:N) {
        g = g + ((( 2 * n + 1) * P[n]) / ((n * n + n)^m    ))
        h = h + (((-2 * n - 1) * P[n]) / ((n * n + n)^(m-1)))
      }
        G[i, j] =  g / 4 / pi  # finalize cell of G-matrix
        H[i, j] = -h / 4 / pi  # finalize cell of H-matrix
      }
   }
   cat("]\n")                  # finalize progress bar
   GH <- list(G = G, H = H)
   return(GH)
}
