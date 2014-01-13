\docType{data}
\name{default.montage.lookup.table}
\alias{default.montage.lookup.table}
\title{Look-up table for 465 standard scalp sites.}
\format{A data frame with 8 columns.}
\source{
  http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/CSDtoolbox.zip
}
\description{
  Look-up table for 330 standard scalp sites based on the
  extended 10-20-system (10/20-, 10/10- or 10/5-system; cf.
  Jurcak et al., 2007), plus additional locations for left
  and right preauricular points (A1/A2; LPA/RPA = T9/T10),
  left and right mastoids (LM/RM = TP9/TP10), and the nose
  tip (Nose).

  In addition, 129 scalp locations of a first-generation
  geodesic sensor net (GSN; Electrical Geodesics, Inc.) are
  also provided. However, in deviation from the original
  129-channel GSN montage, the location for channel #17
  points to the Nose instead of Nz; if appropriate, this
  location must be adjusted in the look-up table by
  replacing the values for channel #17 with the values for
  Nz.
}
\note{
  Unlike the original 129-channel GSN montage, the location
  for channel #17 points to the Nose instead of Nz; if
  appropriate, this location must be adjusted in the
  look-up table by replacing the values for channel #17
  with the values for Nz.
}
\section{Variables}{
  \itemize{

  \item \code{ElectrodeLabel}: Standard electrode labels.

  \item \code{Theta}: rotation of x-axis towards the y-axis
  (in spherical coordinates).

  \item \code{Phi}: rotation of x-y plane towards the
  z-axis (in spherical coordinates).

  \item \code{Radius}: Radius of the sphere.

  \item \code{X}: x-axis position in cartesian coordinates.

  \item \code{Y}: y-axis position in cartesian coordinates.

  \item \code{Z}: z-axis position in cartesian coordinates.

  \item \code{OffSphereSurface}

  }
}
\references{
  Oostenveld, R., Praamstra, P. (2001). The five percent
  electrode system for high-resolution EEG and ERP
  measurements. Clinical Neurophysiology, 112(4), 713-719.

  Jurcak, V., Tsuzuki, D., Dan, I. (2007). 10/20, 10/10,
  and 10/5 systems revisited: their validity as relative
  head-surface-based positioning systems. NeuroImage,
  34(4), 1600-1611.
}
\keyword{datasets}

