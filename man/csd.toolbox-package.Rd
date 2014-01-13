\docType{package}
\name{csd.toolbox-package}
\alias{csd}
\alias{csd.toolbox-package}
\alias{toolbox}
\title{Tools to calculate current source density estimates (CSD) of
  electroencephalographic (EEG) and event-related brain potentials (ERP)
  data}
\description{
  The current source density (CSD) estimates of EEG and ERP
  data can be a useful source of information regarding the
  interpretation of brain activity. The CSD transformation
  provides a reference-free representation of EEG/ERP data
  that reveals primarily local, shallow and
  radially-oriented source activity. A direct comparison
  with EEG or surface potentials can thus help establish
  their source composition: EEG/ERP data will show the
  contribution of shallow and deeper sources as well as
  sources with any orientation, whereas CSD-transformed
  data will isolate local, shallow, radially-oriented
  sources (see Srinivasan, 2005 for example).

  In the context of source localization using equivalent
  current dipole modeling (ECD), this can serve as a
  reasonable check of assumptions for whether (and if so,
  how many) dipolar sources one should posit in their
  source space.

  In the context of concurrent magnetoencephalography (MEG)
  recordings, the CSD representation might be fruitfully
  combined with MEG data, as each provide optimal
  sensitivities to sources oriented orthogonally from each
  other (radial for CSD, tangential for MEG).

  In the context of oscillatory dynamics of brain activity,
  the CSD transformation might provide a more accurate
  record of local vs long-distance cortical synchronization
  and desynchronization, as the effects of volume
  conduction are drastically reduced after the
  transformation.

  Finally, the CSD transformation might be used as an end
  in itself, in cases where it provides a better isolation
  of EEG or evoked activity of interest.

  This package provides tools to calculate the CSD
  transformation of EEG and ERP data. The basic functions
  have been ported from the MATLAB CSD Toolbox created by
  Jürgen Kayser. These functions implement the spherical
  spline algorithm developed by Perrin et al. (1989, 1990).
  The package also adapted from the original MATLAB toolbox
  a look-up table of the location of 330 standard scalp
  sites based on the extended 10-20-system (10/20-, 10/10-
  or 10/5-system; cf. Jurcak et al., 2007), in addition to
  129 scalp locations of a first-generation geodesic sensor
  net.
}
\author{
  Diogo Almeida \email{dalazal@gmail.com}
}
\references{
  Kayser, J., Tenke, C.E. (2006a). Principal components
  analysis of Laplacian waveforms as a generic method for
  identifying ERP generator patterns: I. Evaluation with
  auditory oddball tasks. Clinical Neurophysiology, 117(2),
  348-368.

  Kayser, J., Tenke, C.E. (2006b). Principal components
  analysis of Laplacian waveforms as a generic method for
  identifying ERP generator patterns: II. Adequacy of
  low-density estimates. Clinical Neurophysiology, 117(2),
  369-380.

  Kayser, J. (2009). Current source density (CSD)
  interpolation using spherical splines - CSD Toolbox
  (Version 1.1)
  [http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox].
  New York State Psychiatric Institute: Division of
  Cognitive Neuroscience.

  Perrin, F., Pernier, J., Bertrand, O., Echallier, J.F.
  (1989).  Spherical splines for scalp potential and
  current density mapping.  Electroencephalography and
  Clinical Neurophysiology, 72(2), 184-187.

  Perrin, F., Pernier, J., Bertrand, O., Echallier, J.F.
  (1990).  Corrigenda EEG 02274. Electroencephalography and
  Clinical Neurophysiology, 76, 565.

  Oostenveld, R., Praamstra, P. (2001). The five percent
  electrode system for high-resolution EEG and ERP
  measurements. Clinical Neurophysiology, 112(4), 713-719.

  Jurcak, V., Tsuzuki, D., Dan, I. (2007). 10/20, 10/10,
  and 10/5 systems revisited: their validity as relative
  head-surface-based positioning systems. NeuroImage,
  34(4), 1600-1611.

  Srinivasan, R. (2005). High-Resolution EEG: Theory and
  Practice.  In: Handy, T. (ed.) Event-related potentials:
  A methods handbook, 167-188.
}
\keyword{package}

