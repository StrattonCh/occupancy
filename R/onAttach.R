# USGS require a message on packages that have not been officially released. Adding this to any file in the R folder will cause a message to be displayed when a user opens the library:

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Although this software program has been used by the U.S. Geological Survey (USGS), no warranty, expressed or implied, is made by the USGS or the U.S. Government as to the accuracy and functioning of the program and related program material nor shall the fact of distribution constitute any such warranty, and no responsibility is assumed by the USGS in connection therewith.")
}
