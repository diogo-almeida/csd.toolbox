\name{ExtractMontage}
\alias{ExtractMontage}
\title{Extracts an EEG montage from a standard coordinate look-up table.}
\usage{
  ExtractMontage(user.labels,
    montage.lookup.table = default.montage.lookup.table)
}
\arguments{
  \item{user.labels}{Vector supplying the electrode labels
  used in the user's montage.}

  \item{montage.lookup.table}{Matrix containing the
  information from which the user's montage can be
  extracted. Defaults to matrix
  \code{\link{default.montage.lookup.table}}, which is
  automatically loaded with the package.}
}
\value{
  List containing location and label information about the
  user's montage, with the following fields:

  \describe{

  \item{electrode.labels}{Vector containing the standard
  electrode labels used in the montage.}

  \item{theta}{Vector containing the azimuthal locations of
  the electrodes used in the montage (in spherical
  coordinates).}

  \item{phi}{Vector containing the elevation locations of
  the electrodes used in the montage (in spherical
  coordinates).}

  \item{xy}{Matrix containing the xy cartesian locations of
  the electrodes used in the montage.}

  \item{not.found}{Vector containing electrodes labels
  supplied by the user that have not been found in the look
  up table. If all electrodes were found, it contains NA} }
}
\description{
  Extracts an EEG montage from a standard coordinates
  look-up table based on a vector of labels of standard
  electrode sites. The look up table can be found in the
  \code{default.montage.lookup.table} matrix that is
  automatically loaded with the package. The default look
  up table contains 330 standard 10-20-system and 129
  geodesic (GSN) scalp locations.  Users can supply their
  own look up tables if necessary (see
  \code{\link{default.montage.lookup.table}} for formatting
  details).
}
\examples{
montage <- ExtractMontage(colnames(NR.C66.trr))
str(montage)

electrode.labels <- c("FPz", "Cz", "CPz", "Pz")
montage <- ExtractMontage(electrode.labels)
str(montage)

electrode.labels <- c("FPz", "Cz", "CPz", "Pz", "MR541")
montage <- ExtractMontage(electrode.labels)
str(montage)
}
\seealso{
  \code{\link{default.montage.lookup.table}}
}

