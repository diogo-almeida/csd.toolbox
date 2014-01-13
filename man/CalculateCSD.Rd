\name{CalculateCSD}
\alias{CalculateCSD}
\title{Calculates the Current Source Density estimates of EEG/ERP data.}
\usage{
  CalculateCSD(eeg.data, G, H, head.radius = 1,
    lambda = 10^-5, scalp.potentials = TRUE)
}
\arguments{
  \item{eeg.data}{Matrix containing the data to be
  transformed.}

  \item{G}{The transformation matrix for the interpolation
  of the surface potentials, as calculated by function
  \code{\link{GetGH}}.}

  \item{H}{The transformation matrix for the calculation of
  the Current Source Density (CSD) estimates, as calculated
  by function \code{\link{GetGH}}.}

  \item{head.radius}{Head radius. Defaults to 1 for unit
  sphere [µV/m²]).  Specify a value (in cm) to rescale CSD
  data to smaller units [µV/cm²] (e.g., use 10.0 to scale
  to more realistic head size).}

  \item{lambda}{Smoothing constant lambda for the spherical
  spline interpolation. Defaults to 0.00001.}

  \item{scalp.potentials}{Logical value. \code{TRUE}
  determines the the spherical spline surface potential
  interpolation on top of the calculation of the CSD
  estimates. \code{FALSE} makes the function return only
  the CSD estimates. Defaults to \code{TRUE}.}
}
\value{
  List containing the spherical spline interpolation of the
  CSD and surface potential data (in case the latter was
  requested), with the following fields:

  \describe{

  \item{csd}{Matrix containing the CSD transformed data.}

  \item{surface.potential}{Matrix containing the
  interpolation of the surface potential data, or NA, in
  case this was not requested by the user.}

  }
}
\description{
  Calculates the Current Source Density (CSD) estimates of
  EEG/ERP data by means of the spherical spline
  interpolation algorithm developed by Perrin et al (1989,
  1990).
}
\examples{
\dontrun{
# From the original MATLAB CSD Toolbox tutorial
m.example  <- ExtractMontage(user.labels = colnames(NR.C66.trr))
gh.example <- GetGH(m.example)
csd.data   <- CalculateCSD(eeg.data = t(NR.C66.trr), G = gh.example$G,
  H = gh.example$H
)

# Should reproduce Figures 16 of the original MATLAB CSD TOOLBOX tutorial:
# http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/tutorial.html
matplot(t(csd.data$csd), type = "l", lty = 1, main = "Figure 16",
  ylim = c(-20, 40)
)

# Should reproduce Figures 17 of the original MATLAB CSD TOOLBOX tutorial
matplot(t(csd.data$csd)[, c(14,24)], type = "l", col = c("blue", "green"),
  lty = 1, main = "Figure 17", ylim = c(-10, 35)
 )

# Should reproduce Figures 18 of the original MATLAB CSD TOOLBOX tutorial
matplot(NR.C66.trr, type = "l", lty = 1, main = "Figure 18",
  ylim = c(-10, 20)
)

# Should reproduce Figures 19 of the original MATLAB CSD TOOLBOX tutorial
matplot(NR.C66.trr[, c(14,24)], type = "l", col = c("blue", "green"),
  lty = 1, main = "Figure 19", ylim = c(-10, 20)
)
}
}
\seealso{
  \code{\link{GetGH}}
}

