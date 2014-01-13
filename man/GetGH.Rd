\name{GetGH}
\alias{GetGH}
\title{Computes the transformation matrices G and H for CSD transformation.}
\usage{
  GetGH(M, m = 4)
}
\arguments{
  \item{M}{Montage list (output of ExtractMontage
  function), containing the spherical and cartesian
  coordinates of the electrode montage to be used in the
  CSD transformation.}

  \item{m}{Value of the m-constant to be used in the
  spherical spline interpolation. Lower values allow more
  flexible splines than higher values.  The range of
  possible values goes from 2 to 10, 4 being the default.
  Changing this value changes the properties of the spatial
  filter, so you need to know what you are doing when
  changing this value.}
}
\value{
  List containing the two transformation matrices G and H.
  Matrix G is used in computing the spherical spline
  interpolation of the surface potentials, whereas matrix H
  is used in the spherical spline interpolation that leads
  to the current source density estimate of the data.
}
\description{
  Computes the transformation matrices G and H to be used
  in the CSD transformation.
}
\note{
  This function requires computing Legendre polynomials.
  The original MATLAB code used the function
  \code{legendre}. The R base lacks a similar function, but
  package \code{\link{orthopolynom}} implements an
  equivalent solution, even though the user interface is
  substantially different from the MATLAB original
  function.

  IMPORTANT: The Legendre polynomial results from MATLAB
  and R are not exactly the same all the time. The reasons
  at this point are unclear to me.  However, despite not
  giving exactly the same results, tests based on the
  tutorial data of the original MATLAB CSD Toolbox indicate
  that the MATLAB output and the R output are equal up to
  the 6th and the 10th decimal places for G and H
  respectively, so the differences seem to be negligible.
}
\examples{
montage <- ExtractMontage(colnames(NR.C66.trr))
}
\seealso{
  \code{\link{ExtractMontage}}; \code{\link{CalculateCSD}};
  \code{\link{orthopolynom}}
}

