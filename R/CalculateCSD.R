#' @title Calculates the Current Source Density estimates of EEG/ERP data.
#' 
#' @description Calculates the Current Source Density (CSD) estimates of 
#'   EEG/ERP data by means of the spherical spline interpolation algorithm
#'   developed by Perrin et al (1989, 1990).
#'
#' @param eeg.data Matrix containing the data to be transformed.
#'
#' @param G The transformation matrix for the interpolation of the surface 
#'   potentials, as calculated by function \code{\link{GetGH}}.
#'
#' @param H The transformation matrix for the calculation of the Current Source 
#'   Density (CSD) estimates, as calculated by function \code{\link{GetGH}}.
#'
#' @param head.radius Head radius. Defaults to 1 for unit sphere [µV/m²]).
#'   Specify a value (in cm) to rescale CSD data to smaller units [µV/cm²]
#'   (e.g., use 10.0 to scale to more realistic head size).
#'
#' @param lambda Smoothing constant lambda for the spherical spline 
#'   interpolation. Defaults to 0.00001.
#'
#' @param scalp.potentials Logical value. \code{TRUE} determines the the 
#'   spherical spline surface potential interpolation on top of the calculation 
#'   of the CSD estimates. \code{FALSE} makes the function return only the CSD 
#'   estimates. Defaults to \code{TRUE}.
#'
#' @return List containing the spherical spline interpolation of the CSD and 
#'   surface potential data (in case the latter was requested), with the 
#'   following fields:
#'
#' \describe{
#' 
#' \item{csd}{Matrix containing the CSD transformed data.}
#' 
#' \item{surface.potential}{Matrix containing the interpolation of the surface 
#'   potential data, or NA, in case this was not requested by the user.}
#'
#' }
#'
#' @export
#'
#' @seealso \code{\link{GetGH}}
#'
#' @examples
#' 
#' \dontrun{
#' # From the original MATLAB CSD Toolbox tutorial
#' m.example  <- ExtractMontage(user.labels = colnames(NR.C66.trr))
#' gh.example <- GetGH(m.example)
#' csd.data   <- CalculateCSD(eeg.data = t(NR.C66.trr), G = gh.example$G, 
#'   H = gh.example$H
#' )
#' 
#' # Should reproduce Figures 16 of the original MATLAB CSD TOOLBOX tutorial:
#' # http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/tutorial.html
#' matplot(t(csd.data$csd), type = "l", lty = 1, main = "Figure 16", 
#'   ylim = c(-20, 40)
#' )
#'
#' # Should reproduce Figures 17 of the original MATLAB CSD TOOLBOX tutorial
#' matplot(t(csd.data$csd)[, c(14,24)], type = "l", col = c("blue", "green"), 
#'   lty = 1, main = "Figure 17", ylim = c(-10, 35)
#'  )
#' 
#' # Should reproduce Figures 18 of the original MATLAB CSD TOOLBOX tutorial
#' matplot(NR.C66.trr, type = "l", lty = 1, main = "Figure 18", 
#'   ylim = c(-10, 20)
#' )
#' 
#' # Should reproduce Figures 19 of the original MATLAB CSD TOOLBOX tutorial
#' matplot(NR.C66.trr[, c(14,24)], type = "l", col = c("blue", "green"), 
#'   lty = 1, main = "Figure 19", ylim = c(-10, 20)
#' )
#' }
CalculateCSD <- function(eeg.data, G, H, head.radius = 1, lambda = 10^-5,
                         scalp.potentials = TRUE) {
  nElec <- dim(eeg.data)[1]  # get data matrix dimensions
  nPnts <- dim(eeg.data)[2]
  # Error checking of input matrix dimensions
  if ( 
    (dim(G)[1] != dim(G)[2]) |   # matrix dimension error checking
    (dim(H)[1] != dim(H)[2]) |   # G and H matrix must be nElec-by-nElec
    (dim(G)[1] != nElec) | 
    (dim(H)[1] != nElec)) {
      stop(
        sprintf(
          "G (%d-by-%d) and H (%d-by-%d) matrix, dimensions must match rows (%d) of data matrix", 
          dim(G)[1], dim(G)[2], dim(H)[1], dim(H)[2], nElec
        )
      ) 
  }
  mu <- colMeans(eeg.data)         # get grand mean
  Z  <- sweep(eeg.data, 2, mu)     # compute average reference

  Y <- matrix(numeric(nElec * nPnts), ncol = nPnts) # claim memory for
  X <- Y											# output matrices

  head.radius <- head.radius^2 # scaling variable defaults to 1 : [µV/m²]
                               # rescale data to head sphere if 
                               # head.radius > 1 : [µV/cm²]
  
  for (e in 1:dim(G)[1]) {       # add smoothing constant (lambda) to diagonal 
    G[e, e] <- G[e, e] + lambda  # lambda defaults to 1*(10^-5)
  }

  Gi <- solve(G)               # compute G inverse
  TC <- rowSums(Gi)            # compute sums for each row
  sgi <- sum(TC)               # compute sum total

  for (p in 1:nPnts) {
    Cp <- Gi %*% Z[, p]        # compute preliminary C vector
    c0 <- sum(Cp) / sgi        # common constant across electrodes
    C <- Cp - (c0 * TC)        # compute final C vector
    for (e in 1:nElec) {       # compute all CSDs ...
      X[e, p] <- sum(C * H[e, ]) / head.radius  # ... and scale to head size
    }
    if (scalp.potentials) {    # if requested, compute all SPs
      for (e in 1:nElec) { 
        Y[e, p] <- c0 + sum(C * G[e, ])
      }
    }
  }

  if (scalp.potentials) {
    interpolation <- list(csd = X, scalp.potentials = Y)   
  } else {
    interpolation <- list(csd = X, scalp.potentials = NA)
  }
  return(interpolation)
}
