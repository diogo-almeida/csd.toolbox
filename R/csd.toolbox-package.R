#' @title Tools to calculate current source density estimates (CSD) of 
#'   electroencephalographic (EEG) and event-related brain potentials (ERP)
#'   data
#'   
#' @description The current source density (CSD) estimates of EEG and ERP data 
#'   can be a useful source of information regarding the interpretation of brain 
#'   activity. The CSD transformation provides a reference-free representation 
#'   of EEG/ERP data that reveals primarily local, shallow and radially-oriented
#'   source activity. A direct comparison with EEG or surface potentials can 
#'   thus help establish their source composition: EEG/ERP data will show the 
#'   contribution of shallow and deeper sources as well as sources with any 
#'   orientation, whereas CSD-transformed data will isolate local, shallow, 
#'   radially-oriented sources (see Srinivasan, 2005 for example).
#'
#'   In the context of source localization using equivalent current dipole 
#'   modeling (ECD), this can serve as a reasonable check of assumptions for 
#'   whether (and if so, how many) dipolar sources one should posit in their 
#'   source space.
#'
#'   In the context of concurrent magnetoencephalography (MEG) recordings, the 
#'   CSD representation might be fruitfully combined with MEG data, as each 
#'   provide optimal sensitivities to sources oriented orthogonally from each 
#'   other (radial for CSD, tangential for MEG).
#'
#'   In the context of oscillatory dynamics of brain activity, the CSD 
#'   transformation might provide a more accurate record of local vs 
#'   long-distance cortical synchronization and desynchronization, as the 
#'   effects of volume conduction are drastically reduced after the 
#'   transformation.
#'   
#'   Finally, the CSD transformation might be used as an end in itself, in cases
#'   where it provides a better isolation of EEG or evoked activity of interest.
#'
#'   This package provides tools to calculate the CSD transformation of EEG and 
#'   ERP data. The basic functions have been ported from the MATLAB CSD Toolbox
#'   created by Jürgen Kayser. These functions implement the spherical spline 
#'   algorithm developed by Perrin et al. (1989, 1990). The package also 
#'   adapted from the original MATLAB toolbox a look-up table of the location 
#'   of 330 standard scalp sites based on the extended 10-20-system (10/20-, 
#'   10/10- or 10/5-system; cf. Jurcak et al., 2007), in addition to 129 scalp 
#'   locations of a first-generation geodesic sensor net.
#'
#' @references Kayser, J., Tenke, C.E. (2006a). Principal components analysis of
#'   Laplacian waveforms as a generic method for identifying ERP generator 
#'   patterns: I. Evaluation with auditory oddball tasks. Clinical 
#'   Neurophysiology, 117(2), 348-368.
#' 
#' @references Kayser, J., Tenke, C.E. (2006b). Principal components analysis of
#'   Laplacian waveforms as a generic method for identifying ERP generator 
#'   patterns: II. Adequacy of low-density estimates. Clinical Neurophysiology, 
#'   117(2), 369-380.
#'
#' @references Kayser, J. (2009). Current source density (CSD) interpolation 
#'   using spherical splines - CSD Toolbox (Version 1.1) 
#'   [http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox]. New York 
#'   State Psychiatric Institute: Division of Cognitive Neuroscience.
#'
#' @references Perrin, F., Pernier, J., Bertrand, O., Echallier, J.F. (1989). 
#'   Spherical splines for scalp potential and current density mapping. 
#'   Electroencephalography and Clinical Neurophysiology, 72(2), 184-187.
#' 
#' @references Perrin, F., Pernier, J., Bertrand, O., Echallier, J.F. (1990). 
#'   Corrigenda EEG 02274. Electroencephalography and Clinical Neurophysiology, 
#'   76, 565.
#'
#' @references Oostenveld, R., Praamstra, P. (2001). The five percent electrode 
#'   system for high-resolution EEG and ERP measurements. Clinical 
#'   Neurophysiology, 112(4), 713-719.
#'
#' @references Jurcak, V., Tsuzuki, D., Dan, I. (2007). 10/20, 10/10, and 10/5 
#'   systems revisited: their validity as relative head-surface-based 
#'   positioning systems. NeuroImage, 34(4), 1600-1611.
#'
#' @references Srinivasan, R. (2005). High-Resolution EEG: Theory and Practice. 
#'   In: Handy, T. (ed.) Event-related potentials: A methods handbook, 167-188.
#'
#' @name csd.toolbox-package
#' @aliases csd toolbox
#' @docType package
#' @author Diogo Almeida \email{dalazal@@gmail.com}
#' @keywords package
NULL
#' @title Look-up table for 465 standard scalp sites.
#'
#' @description Look-up table for 330 standard scalp sites based on the 
#'   extended 10-20-system (10/20-, 10/10- or 10/5-system; cf. Jurcak et al., 
#'   2007), plus additional locations for left and right preauricular points
#'   (A1/A2; LPA/RPA = T9/T10), left and right mastoids (LM/RM = TP9/TP10), and 
#'   the nose tip (Nose).
#'
#'   In addition, 129 scalp locations of a first-generation geodesic sensor net
#'   (GSN; Electrical Geodesics, Inc.) are also provided. However, in deviation
#'   from the original 129-channel GSN montage, the location for channel 
#'   #17 points to the Nose instead of Nz; if appropriate, this location must be
#'   adjusted in the look-up table by replacing the values for channel #17 with 
#'   the values for Nz.
#' 
#' @name default.montage.lookup.table
#' @docType data
#' @usage default.montage.lookup.table
#' @format A data frame with 8 columns.
#' @section Variables:
#'
#' \itemize{
#' 
#' \item \code{ElectrodeLabel}: Standard electrode labels.
#' 
#' \item \code{Theta}: rotation of x-axis towards the y-axis (in spherical 
#'   coordinates).
#' 
#' \item \code{Phi}: rotation of x-y plane towards the z-axis (in spherical 
#'   coordinates).
#' 
#' \item \code{Radius}: Radius of the sphere.
#' 
#' \item \code{X}: x-axis position in cartesian coordinates.
#' 
#' \item \code{Y}: y-axis position in cartesian coordinates.
#' 
#' \item \code{Z}: z-axis position in cartesian coordinates.
#' 
#' \item \code{OffSphereSurface}
#'
#' }
#'
#' @note Unlike the original 129-channel GSN montage, the location for channel 
#'   #17 points to the Nose instead of Nz; if appropriate, this location must be
#'   adjusted in the look-up table by replacing the values for channel #17 with 
#'   the values for Nz.
#'
#' @references Oostenveld, R., Praamstra, P. (2001). The five percent electrode 
#'   system for high-resolution EEG and ERP measurements. Clinical 
#'   Neurophysiology, 112(4), 713-719.
#'
#' @references Jurcak, V., Tsuzuki, D., Dan, I. (2007). 10/20, 10/10, and 10/5 
#'   systems revisited: their validity as relative head-surface-based 
#'   positioning systems. NeuroImage, 34(4), 1600-1611.
#'
#' @source http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/CSDtoolbox.zip
NULL
#' @title Grand-average ERP data from 66 participants in auditory oddball task.
#'
#' @description Grand-average ERP data derived from 66 healthy adults that
#'   participated in an auditory oddball task (Kayser & Tenke, 2006a). The data 
#'   refers to their complex tones condition, where participants were instructed 
#'   to press the right button when rare stimuli were perceived. The data was 
#'   acquired with a 31 electrode montage, at a sampling frequency of 200 Hz. 
#'   This dataset has 200 samples, ranging from time 0 to 995 ms post-stimulus 
#'   onset, and has been referenced to the tip of the nose.
#' 
#' @name NR.C66.trr
#' @docType data
#' @usage NR.C66.trr
#' @format A 200 x 31 matrix. Rows are time samples (in increments of 5 ms). 
#'   Columns are electrode sites.
#'
#' @references Kayser, J., Tenke, C.E. (2006a). Principal components analysis of
#'   Laplacian waveforms as a generic method for identifying ERP generator 
#'   patterns: I. Evaluation with auditory oddball tasks. Clinical 
#'   Neurophysiology, 117(2), 348-368.
#'
#' @source http://psychophysiology.cpmc.columbia.edu/software/CSDtoolbox/CSDtoolbox.zip
NULL