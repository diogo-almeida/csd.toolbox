#' @title Extracts an EEG montage from a standard coordinate look-up table.
#' 
#' @description Extracts an EEG montage from a standard coordinates look-up 
#'   table based on a vector of labels of standard electrode sites. The look up
#'   table can be found in the \code{default.montage.lookup.table} matrix that
#'   is automatically loaded with the package. The default look up table 
#'   contains 330 standard 10-20-system and 129 geodesic (GSN) scalp locations.
#'   Users can supply their own look up tables if necessary (see 
#'   \code{\link{default.montage.lookup.table}} for formatting details).
#'
#' @param user.labels Vector supplying the electrode labels used in the user's
#'   montage.
#'
#' @param montage.lookup.table Matrix containing the information from which the
#'   user's montage can be extracted. Defaults to matrix 
#'   \code{\link{default.montage.lookup.table}}, which is automatically loaded
#'   with the package.
#'
#' @return List containing location and label information about the user's 
#'   montage, with the following fields:
#'
#' \describe{
#' 
#' \item{electrode.labels}{Vector containing the standard electrode labels 
#'   used in the montage.}
#' 
#' \item{theta}{Vector containing the azimuthal locations of the electrodes used
#'   in the montage (in spherical coordinates).}
#' 
#' \item{phi}{Vector containing the elevation locations of the electrodes used
#'   in the montage (in spherical coordinates).}
#' 
#' \item{xy}{Matrix containing the xy cartesian locations of the electrodes 
#'   used in the montage.}
#'
#' \item{not.found}{Vector containing electrodes labels supplied by the user
#'   that have not been found in the look up table. If all electrodes were 
#'   found, it contains NA}
#' }
#'
#' @export
#'
#' @seealso \code{\link{default.montage.lookup.table}}
#'
#' @examples
#'
#' montage <- ExtractMontage(colnames(NR.C66.trr))
#' str(montage)
#'
#' electrode.labels <- c("FPz", "Cz", "CPz", "Pz")
#' montage <- ExtractMontage(electrode.labels)
#' str(montage)
#'
#' electrode.labels <- c("FPz", "Cz", "CPz", "Pz", "MR541")
#' montage <- ExtractMontage(electrode.labels)
#' str(montage)
ExtractMontage <-function(user.labels,
                          montage.lookup.table = default.montage.lookup.table) {
                          
  # is electrode.labels a list?
  if (is.list(user.labels) == T) {
    user.labels <- as.vector(unlist(user.labels))
    cat("Converting electrode labels from list to vector format... DONE!\n\n")
  }

  # Get the look up table for electrode names and their Theta and Phi locations
  electrode.labels  <- as.vector(montage.lookup.table$ElectrodeLabel)
  theta <- montage.lookup.table$Theta
  phi   <- montage.lookup.table$Phi

  # Find the indices of the user labels found in the look up table as well as
  # the indices (in user.labels) of the electrodes that might not have been 
  # found.
  idx.elec.found     <- toupper(electrode.labels) %in% toupper(user.labels)
  idx.elec.not.found <- !(toupper(user.labels) %in% toupper(electrode.labels))

  # Get the labels of the electrodes supplied by the user from the look up table
  electrodes.found      <- electrode.labels[idx.elec.found]
  electrodes.not.found  <- user.labels[idx.elec.not.found]

  theta                <- theta[idx.elec.found]
  phi                  <- phi[idx.elec.found]
  correct.user.labels  <- electrode.labels[idx.elec.found]

  reorder.vector      <- order(unlist(sapply(toupper(electrodes.found), 
    function(x) {which(toupper(x) == toupper(user.labels))}))
  )
  electrodes.found    <- electrodes.found[reorder.vector]
  theta               <- theta[reorder.vector]
  phi                 <- phi[reorder.vector]

  idx.user.labels.found <- which(toupper(user.labels) %in% electrodes.found)
  correct.user.labels <- electrode.labels[idx.elec.found]

  
  phiT         <- 90 - phi                    # calculate phi from top of sphere
  theta2       <- (2 * pi * theta) / 360      # convert degrees to radians
  phi2         <- (2 * pi * phiT) / 360
  xycoord      <- pol2cart(theta2, phi2)      # get plane coordinates
  xy           <- cbind(xycoord$x, xycoord$y)
  colnames(xy) <- c("x", "y")
  xy           <- xy / max(xy, na.rm = T)     # set maximum to unit length
  xy           <- xy / 2 + 0.5                # adjust to range 0-1

  if (length(user.labels) != length(electrodes.found)){

    cat("\n=========================================\n",
        "ATTENTION:\n",
        "Only", length(correct.user.labels), "out of", length(user.labels),
        "EEG channel labels you provided were identified!\n",
        "The following channels were not found in the look up table:\n",
        electrodes.not.found,
       "\n=========================================\n"
    )
    extracted.montage <- list(electrode.labels = electrodes.found, theta = theta, 
      phi = phi, xy = xy, not.found = electrodes.not.found
    )
  } else {
    extracted.montage <- list(electrode.labels = electrodes.found, theta = theta, 
      phi = phi, xy = xy, not.found = NA
    )
  }
  return(extracted.montage)
}
