#' Print method for `occ_mod` class
#'
#' @param x An object of class [occ_mod]
#' @param ... Other arguments passed to or from other methods
#'
#' @export print.occ_mod
#' @export

print.occ_mod <- function(x, ...) {
  cat("Single species, single season site occupancy model fit using NIMBLE.\n")
  cat("\nThe following model was fit:\n")
  cat(attr(x, "code"))
  cat("\n\nThe following MCMC algorithm was used:\n")
  cat(paste0(attr(x, "mcmcinfo"), sep = "\n"))
}
