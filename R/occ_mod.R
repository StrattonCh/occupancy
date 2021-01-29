#'@title Fit the single species, single season site occupancy model.
#'
#'@description This function fits the single species, single
#'  season site occupancy model first developed by
#'  \href{https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282002%29083%5B2248%3AESORWD%5D2.0.CO%3B2}{MacKenzie
#'   et al. (2002)}.
#'
#'@details This function fits the single season, single species site occupancy
#'  model using the logit link function. The \code{data} should contain columns
#'  named \code{site}, \code{visit}, and \code{y}. See examples.
#'
#'@param occupancy model declaration for the occupancy portion of the model
#'  using standard linear model syntax
#'@param detection model declaration for the detection portion of the model
#'  using standard linear model syntax
#'@param data data from which to create the detection and covariate matrices.
#'  Each row should represent a single visit to a site. See details and examples
#'  for proper formatting
#'@param niter number of MCMC iterations
#'@param seed optional seed
#'@param save_model logical; should a text file containing the model be
#'  exported?
#'@param model_name character string defining the name of the text file
#'  describing the model if \code{save_model = TRUE}
#' @param nchains number of MCMC chains
#'
#'@example examples/occ_mod_ex.R
#'
#'@return an object of class \code{list} containing the following:\cr
#'  * \code{samples} object of class \code{list} of length \code{nchains}, each
#'  containing a \code{matrix} posterior samples
#'  * \code{loglik} object of class \code{list} of length code{nchains}, each
#'  containing a \code{matrix} of samples of the log posterior likelihood
#'
#'@export
#'
#'@md

occ_mod <- function(occupancy, detection, data, niter = 1000, nchains = 3, seed = NULL,
                    save_model = FALSE, model_name = paste0("occ_model_", Sys.Date())){
  # fix unbound variables issues
  site <- visit <- y <- `1` <- `.` <- NULL

  # set seed
  if(!is.null(seed)) set.seed(seed)

  # convenience
  `%notin%` <- Negate("%in%")

  # check names on data
  if("site" %notin% names(data)) stop("Data must contain column named 'site'")
  if("visit" %notin% names(data)) stop("Data must contain column named 'visit'")
  if("y" %notin% names(data)) stop("Data must contain column named 'y'")

  # model formulas
  occupancy_mod <- stats::as.formula(occupancy)
  detection_mod <- stats::as.formula(detection)

  # model matrices
  ## occupancy
  occupancy_mod_matrix <- stats::model.matrix(object = occupancy_mod, data = data)
  occ_num_params <- ncol(occupancy_mod_matrix)
  occupancy_mod_matrix <- cbind(
    occupancy_mod_matrix, matrix(data$site, ncol = 1, dimnames = list(NULL, "site_"))
  ) %>% unique(.)
  if(nrow(occupancy_mod_matrix) != length(unique(data$site))) stop("There are more unique site-level covariate values than there are sites.")
  if(occ_num_params == 1){
    occupancy_mod_matrix <- as.matrix(occupancy_mod_matrix[,-which(colnames(occupancy_mod_matrix) == "site_")], ncol = 1)
    colnames(occupancy_mod_matrix) <- "(Intercept)"
  } else{
    occupancy_mod_matrix <- occupancy_mod_matrix[,-which(colnames(occupancy_mod_matrix) == "site_")]
  }

  ## detection
  detection_mod_matrix <- stats::model.matrix(object = detection_mod, data = data)

  # prep data for NIMBLE
  nimble_data <- list()
  ## response
  nimble_data$Y <- data %>%
    dplyr::select(site, visit, y) %>%
    dplyr::group_by(site) %>%
    tidyr::pivot_wider(names_from = "visit", values_from = "y") %>%
    dplyr::ungroup() %>%
    dplyr::select(`1`:ncol(.)) %>%
    as.matrix()

  ## occupancy level
  occupancy_covariates <- list()
  for(i in 1:ncol(occupancy_mod_matrix)){
    occupancy_covariates[[colnames(occupancy_mod_matrix)[i]]] <- unname(occupancy_mod_matrix[,i])
  }
  nimble_data$occupancy_covariates <- occupancy_covariates

  ## detection level
  detection_covariates <- list()
  for(i in 1:ncol(detection_mod_matrix)){
    detection_covariates[[colnames(detection_mod_matrix)[i]]] <- data %>%
      dplyr::select(site, visit) %>%
      dplyr::mutate(tmp = detection_mod_matrix[,i]) %>%
      dplyr::group_by(site) %>%
      tidyr::pivot_wider(names_from = "visit", values_from = "tmp") %>%
      dplyr::ungroup() %>%
      dplyr::select(`1`:ncol(.)) %>%
      as.matrix()
  }
  nimble_data$detection_covariates <- detection_covariates

  ## nvisits
  nimble_data$nvisits <- apply(nimble_data$Y, 1, function(x) sum(!is.na(x)))

  # NIMBLE code
  ## priors
  priors <- paste0(
    "# priors\n",
    "for(i in 1:", ncol(occupancy_mod_matrix), "){beta_psi[i] ~ dnorm(0, 1/2)}",
    "\n",
    "for(i in 1:", ncol(detection_mod_matrix), "){beta_p[i] ~ dnorm(0, 1/2)}",
    "\n"
  )

  ## likelihood
  likelihood_a <- paste0(
    "# likelihood \n",
    "# loop through sites \n",
    "for(i in 1:M){ \n",
    " # calculate psi \n")
  if(ncol(occupancy_mod_matrix) > 1){
    likelihood_b <- paste0(
      " logit(psi[i]) <- beta_psi[1] +\n",
      paste0(sapply(2:ncol(occupancy_mod_matrix), function(x){
        paste0("  beta_psi[", x, "] * ", colnames(occupancy_mod_matrix)[x], "[i]")
      }), collapse = " + \n")
    )
  } else{
    likelihood_b <- paste0(
      " logit(psi[i]) <- beta_psi[1]"
    )
  }
  likelihood_c <- paste0(
    "\n\n",
    " # latent occupancy state \n",
    " z[i] ~ dbern(psi[i]) \n\n",
    " # loop through visits \n",
    " for(j in 1:nvisits[i]){\n",
    "  # calculate p\n"
  )
  if(ncol(detection_mod_matrix) > 1){
    likelihood_d <- paste0(
      "  logit(p[i,j]) <- beta_p[1] +\n",
      paste0(sapply(2:ncol(detection_mod_matrix), function(x){
        paste0("   beta_p[", x, "] * ", colnames(detection_mod_matrix)[x], "[i,j]")
      }), collapse = " + \n")
    )
  } else{
    likelihood_d <- paste0(
      "  logit(p[i,j]) <- beta_p[1]"
    )
  }
  likelihood_e <- paste0(
    "\n\n",
    "  # response model\n",
    "  Y[i,j] ~ dbern(z[i] * p[i,j])",
    "\n }",
    "\n}"
  )
  likelihood <- paste0(likelihood_a, likelihood_b, likelihood_c, likelihood_d, likelihood_e)

  code <- paste0(priors, "\n", likelihood)

  # initialization
  mod_inits <- function(x){
    tmp <- x
    list(
      z = apply(nimble_data$Y, 1, function(x) ifelse(sum(x, na.rm = T) == 0, 0, 1)),
      beta_psi = stats::rnorm(ncol(occupancy_mod_matrix), 0, sqrt(2)),
      beta_p = stats::rnorm(ncol(detection_mod_matrix), 0, sqrt(2))
    )
  }

  # model
  data_ <- c(nimble_data$occupancy_covariates, nimble_data$detection_covariates)
  data_ <- data_[-c(which(names(data_) == "(Intercept)"))]
  data_[["M"]] <- nrow(nimble_data$Y)
  data_[["nvisits"]] <- nimble_data$nvisits
  data_[["Y"]] <- nimble_data$Y

  # write model
  if(save_model) message("Writing NIMBLE model to current working directory:")
  cat(code)

  fileConn <- file(paste0(model_name, ".txt"))
  writeLines(code, fileConn)
  close(fileConn)

  message("\nBuilding R model")
  model <- nimble::readBUGSmodel(
    model = paste0(model_name, ".txt"),
    data = data_,
    inits = mod_inits(1)
  )
  if(!save_model) file.remove(paste0(model_name, ".txt"))

  # fit model
  message("\nCompiling R model to C")
  model_c <- nimble::compileNimble(model)

  message("\nBuilding R MCMC object")
  mcmcinfo <- capture.output(mcmc <- nimble::buildMCMC(model_c, monitors = c("z", "beta_psi", "beta_p")))

  message("\nCompiling R MCMC object to C")
  mcmc_c <- nimble::compileNimble(mcmc)

  message("\nRunning MCMC")
  inits <- lapply(1:nchains, function(x) mod_inits(x))
  samples <- nimble::runMCMC(mcmc_c, niter = niter, inits = inits, nchains = nchains)

  # compute posterior likelihood
  message("\nComputing posterior likelihood")
  loglik <- list()
  for(i in 1:nchains){
    samples_ <- samples[[i]]
    beta_psi <- samples_[,which(colnames(samples_) %in% paste0('beta_psi[', 1:ncol(occupancy_mod_matrix), ']'))]
    beta_p <- samples_[,which(colnames(samples_) %in% paste0('beta_p[', 1:ncol(detection_mod_matrix), ']'))]
    z <- samples_[,which(colnames(samples_) %in% paste0('z[', 1:nrow(data_[["Y"]]), ']'))]

    z_ <- z[,rep(1:ncol(z),data_[["nvisits"]])]
    y_ <- matrix(rep(data[,"y"],  nrow(samples_)), nrow = nrow(samples_), byrow = TRUE)

    psi_tmp <- exp(beta_psi %*% t(occupancy_mod_matrix))
    psi <- psi_tmp / (1 + psi_tmp)
    p_tmp <- exp(beta_p %*% t(detection_mod_matrix))
    p <- p_tmp / (1 + p_tmp)

    # log pointwise predictive density
    site <- psi^z * (1-psi)^(1-z)

    visit <- (z_*p)^y_ * (1-z_*p)^(1-y_)
    ndx_upr <- cumsum(data_[["nvisits"]])
    ndx_lwr <- ndx_upr - data_[["nvisits"]] + 1
    tmp <- sapply(1:ncol(site), function(x){ apply(visit[,ndx_lwr[x]:ndx_upr[x]], 1, prod)})

    loglik[[i]] <- log(site * tmp)
    colnames(loglik[[i]]) <- paste0("site", 1:nrow(data_[["Y"]]))
  }
  names(loglik) <- paste0("chain", 1:nchains)

  # return
  out <- list(samples = samples, loglik = loglik)

  # attributes
  class(out) <- c("occ_mod", "list")
  attr(out, "code") <- code
  attr(out, "mcmcinfo") <- mcmcinfo
  attr(out, "occupancy_call") <- occupancy
  attr(out, "detection_call") <- detection
  return(out)
}

#' Print method for \code{occ_mod} class
#'
#' @param x An object of class \link{occ_mod}
#' @param ... Other arguments passed to or from other methods
#'
#' @export
print.occ_mod <- function(x, ...) {
  cat("Single species, single season site occupancy model fit using NIMBLE.\n")
  cat("\nThe following model was fit:\n")
  cat(attr(x, "code"))
  cat("\n\nThe following MCMC algorithm was used:\n")
  cat(paste0(attr(x, "mcmcinfo"), sep = "\n"))
}

#' Summary method for \code{occ_mod} class
#'
#' @param x An object of class \link{occ_mod}
#' @param burnin number of iterations to discard as burnin
#' @param waic_type Type of WAIC to calculate (either 1 or 2)
#' @param digits number of digits to print in summary
#' @param ... Other arguments passed to or from other methods
#'
#' @export
summary.occ_mod <- function(x, burnin = nrow(x[["samples"]][[1]])/2, waic_type = 2,
                            digits = max(3L, getOption("digits") - 3L), ...){
  # housekeeping
  niter <- nrow(x[["samples"]][[1]])
  nchains <- length(x[["samples"]])

  # grab samples
  samples <- x[["samples"]]
  loglik <- x[["loglik"]]

  # psrf and mpsrf
  reg_coefs <- lapply(samples, function(x) x[,c(which(grepl("_psi[[]", colnames(x))), which(grepl("_p[[]", colnames(x))))])
  gelman <- coda::gelman.diag(coda::as.mcmc.list(lapply(reg_coefs, function(x) coda::as.mcmc(x))))

  # combine chains
  samples_burnin <- do.call("rbind", lapply(samples, function(x) x[(burnin+1):nrow(x),]))
  loglik_burnin <- do.call("rbind", lapply(loglik, function(x) x[(burnin+1):nrow(x),]))

  # create call
  call <- paste0(
    "occ_mod(occupancy = ",
    paste0(as.character(attr(x, "occupancy_call")), collapse = ""),
    ", detection = ",
    paste0(as.character(attr(x, "detection_call")), collapse = ""),
    ")"
  )

  # occupancy coefficients
  occupancy_coef <- samples_burnin[, grepl("_psi[[]", colnames(samples_burnin))]
  occ_coef <- cbind(
    estimate = colMeans(occupancy_coef),
    sd = apply(occupancy_coef, 2, sd),
    `0.025` = apply(occupancy_coef, 2, quantile, prob = .025),
    `0.5` = apply(occupancy_coef, 2, quantile, prob = .5),
    `0.975` = apply(occupancy_coef, 2, quantile, prob = .975)
  )

  # detection coefficients
  detection_coef <- samples_burnin[, grepl("_p[[]", colnames(samples_burnin))]
  detect_coef <- cbind(
    estimate = colMeans(detection_coef),
    sd = apply(detection_coef, 2, sd),
    `0.025` = apply(detection_coef, 2, quantile, prob = .025),
    `0.5` = apply(detection_coef, 2, quantile, prob = .5),
    `0.975` = apply(detection_coef, 2, quantile, prob = .975)
  )

  # waic
  lppd <- sum(log(colMeans(exp(loglik_burnin))))
  pwaic1 <- 2 * sum(log(colMeans(exp(loglik_burnin))) - colMeans(loglik_burnin))
  pwaic2 <- sum(apply(loglik_burnin, 2, var))
  waic1 <- -2 * (lppd - pwaic1)
  waic2 <- -2 * (lppd - pwaic2)

  out <- list(
    call = call,
    occ_coef = occ_coef,
    detect_coef = detect_coef,
    waic1 = waic1,
    waic2 = waic2,
    gelman.diag = gelman
  )
  class(out) <- "summary.occ_mod"
  attr(out, "waic_type") <- waic_type
  attr(out, "digits") <- digits
  return(out)
}

#' Print method for \code{summary.occ_mod} class
#'
#' @param x An object of class \link{summary.occ_mod}
#' @param ... Other arguments passed to or from other methods
#'
#' @export
print.summary.occ_mod <- function(x, ...){
  # call
  cat(
    paste0(
      "Call:\n",
      x[["call"]], "\n"
    )
  )

  # convergence diagnostics
  cat(
    paste0(
      "\nMultivariate potential scale reduction factor:\n ",
      round(x[["gelman.diag"]][['mpsrf']], digits = attr(x, "digits")), "\n"
    )
  )
  cat(
    paste0(
      "\nPotential scale reduction factors:\n"
    )
  )
  print(noquote(t(apply(x[["gelman.diag"]][["psrf"]], 1, formatC, digits = attr(x, "digits")))))

  # occupancy coeff
  cat(
    paste0(
      "\nOccupancy coefficients:\n"
    )
  )
  print(noquote(t(apply(x[["occ_coef"]], 1, formatC, digits = attr(x, "digits")))))

  # detection coef
  cat(
    paste0(
      "\nDetection coefficients:\n"
    )
  )
  print(noquote(t(apply(x[["detect_coef"]], 1, formatC, digits = attr(x, "digits")))))

  # waic
  waic <- ifelse(attr(x, "waic_type") == 1, x[["waic1"]], x[["waic2"]])
  cat(
    paste0(
      "\nType ", attr(x, "waic_type"), " waic: ", round(waic, digits = attr(x, "digits")), "\n"
    )
  )
}



