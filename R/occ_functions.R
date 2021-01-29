#'@title Simulate data from the single species, single season site occupancy model.
#'
#'@description This function simulates data from the single species, single
#'  season occupancy model first developed by
#'  [MacKenzie
#'   et al. (2002)](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282002%29083%5B2248%3AESORWD%5D2.0.CO%3B2).
#'
#' @details This function simulates data from the vanilla single season, single
#'  species occupancy model using the logit link function. If `rand_visits
#'  = TRUE`, each site is visited a random number of times between two and
#'  `max_j`. Covariates are drawn from the uniform(0, 1.5) distribution so
#'  that the effect of the direction of each regression coefficient is
#'  intuitive. Note that if covariates are not desired, `beta_psi` and
#'  `beta_p` can be set to intercepts that generate the desired derived
#'  probabilities.
#'
#'@param M number of sites
#'@param max_j maximum number of visits to each site. If `rand_visits =
#'  FALSE`, this value is the number of visits to each site. If
#'  `rand_visits = TRUE`, each site is visited a random
#'  (`sample(2:max_j, size = 1)`) number of times.
#'@param beta_psi vector of regression coefficients used to generate psi
#'@param beta_p vector of regression coefficients used to generate p
#'@param seed optional seed for reproducibility
#'@param rand_visits logical; should each site be visited a random number of
#'  times? See details.
#'
#'@example examples/sim_occ_ex.R
#'
#'@return object of class `list` containing the following elements: \cr
#'  * `beta_psi` vector of regression coefficients used to generate psi
#'  * `beta_p` vector of regression coefficients used to generate p
#'  * `psi_cov` matrix of site level covariates
#'  * `p_cov` array of detection level covariates; each slice represents a
#'  single covariate
#'  * `psi` vector of derived site level occupancy probabilities
#'  * `p` matrix of derived visit level detection probabilities
#'  * `z` vector of latent occupancy states for each site
#'  * `Y` matrix of observed Bernoulli responses
#'  * `n_visits` vector of number of visits to each site
#'  * `data` a data frame containing all information necessary to fit the
#'  model
#'
#'@importFrom magrittr %>%
#'@export
#'
#'@md

sim_occ <- function(M = 20, max_j = 10, beta_psi = c(0, 1), beta_p = c(0, 1),
                    seed = NULL, rand_visits = TRUE) {
  # create out vector
  out <- list()

  # optional seed
  if (!is.null(seed)) set.seed(seed)

  # convenience
  inv.logit <- function(x) exp(x) / (1 + exp(x))

  # double check betas
  beta_psi <- as.matrix(beta_psi, ncol = 1)
  beta_p <- as.matrix(beta_p, ncol = 1)

  # dimension of regression coefficients
  p_beta_psi <- nrow(beta_psi)
  p_beta_p <- nrow(beta_p)

  # generate covariates
  ## psi
  if(p_beta_psi == 1){
    psi_cov <- matrix(rep(1, M), ncol = 1)
  } else{
    psi_cov <- cbind(
      rep(1, M),
      matrix(
        stats::runif((p_beta_psi-1) * M, 0, 1.5), ncol = p_beta_psi-1
      )
    )
  }

  # column names
  if (p_beta_psi == 1) {
    colnames(psi_cov) <- "psi_int"
  } else {
    colnames(psi_cov) <- c("psi_int", paste0("psi_cov", 1:(ncol(psi_cov) - 1)))
  }

  # generate psi
  psi <- inv.logit(psi_cov %*% beta_psi)

  ## p
  p_cov <- array(0, dim = c(M, max_j, p_beta_p))
  p_linpred <- matrix(0, M, max_j)
  for (k in 1:p_beta_p) {
    if (k == 1) {
      p_cov[, , k] <- 1
    } else {
      p_cov[, , k] <- stats::runif(n = M * max_j, 0, 1.5)
    }

    p_linpred <- p_linpred + beta_p[k, ] * p_cov[, , k]
  }

  # names
  if (p_beta_p == 1) {
    dimnames(p_cov)[[3]] <- "p_int"
  } else {
    dimnames(p_cov)[[3]] <- c("p_int", paste0("p_cov", 1:(dim(p_cov)[3] - 1)))
  }

  # generate p
  p <- inv.logit(p_linpred)

  # generate latent occupancy
  z <- apply(psi, 1, function(x) stats::rbinom(1, 1, x))

  # generate responses
  Y <- matrix(0, M, max_j)
  for (i in 1:M) {
    for (j in 1:max_j) {
      Y[i, j] <- z[i] * stats::rbinom(1, 1, prob = p[i, j])
    }
  }

  # add in NAs if rand_visits
  if (rand_visits) {
    n_visits <- sample(c(2:(max_j)), size = M, replace = TRUE)
    for (i in 1:M) {
      if (n_visits[i] == max_j) {
        NULL
      } else {
        # response
        Y[i, (n_visits[i] + 1):max_j] <- NA

        # covariates
        p_cov[i, (n_visits[i] + 1):max_j, ] <- NA

        # derived parameters
        p[i, (n_visits[i] + 1):max_j] <- NA
      }
    }

    out$n_visits <- apply(Y, 1, function(x) length(which(!is.na(x))))
  }
  if(!rand_visits) out$n_visits <- rep(max_j, M)

  # return
  out$beta_psi <- beta_psi
  out$beta_p <- beta_p
  out$psi_cov <- psi_cov
  out$p_cov <- p_cov
  out$psi <- psi
  out$p <- p
  out$z <- z
  out$Y <- Y

  # create data compatible with occ_mod
  data <- data.frame(
    site = rep(1:M, out$n_visits),
    visit = unlist(c(sapply(out$n_visits, function(x) 1:x))),
    y = stats::na.omit(c(t(out$Y)))
  )

  # site covariates
  tmp <- out$psi_cov
  tmp[,1] <- 1:nrow(tmp);colnames(tmp)[1] <- "site"
  tmp <- as.data.frame(tmp)
  data <- dplyr::left_join(data, tmp, by = "site")

  # detection covariates
  tmp <- as.data.frame(apply(out$p_cov, 3, function(x) stats::na.omit(c(t(x)))))
  tmp$site <- data$site
  tmp$visit <- data$visit
  tmp <- tmp[,-which(names(tmp) == "p_int")]
  data <- dplyr::left_join(data, tmp, by = c("site", "visit"))

  # export
  out$data <- data

  return(out)
}

#'@title Fit the single species, single season site occupancy model.
#'
#'@description This function fits the single species, single
#'  season site occupancy model first developed by
#'  [MacKenzie
#'   et al. (2002)](https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282002%29083%5B2248%3AESORWD%5D2.0.CO%3B2).
#'
#'@details This function fits the single season, single species site occupancy
#'  model using the logit link function. The `data` should contain columns
#'  named `site`, `visit`, and `y`. See examples.
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
#'  describing the model if `save_model = TRUE`
#'@param nchains number of MCMC chains
#'@param beta_prior character string defining prior distribution for regression
#'  coefficients at the occupancy and detection levels. Priors should be
#'  specified using distributions available in NIMBLE. See
#'  [Available
#'   distributons in NIMBLE](https://r-nimble.org/html_manual/cha-writing-models.html#subsec:dists-and-functions). et al. (2002)}.
#'
#'@example examples/occ_mod_ex.R
#'
#'@return an object of class `list` containing the following:\cr
#'  * `samples` object of class `list` of length `nchains`, each
#'  containing a `matrix` posterior samples
#'  * `loglik` object of class `list` of length code{nchains}, each
#'  containing a `matrix` of samples of the log posterior likelihood
#'
#'@export
#'
#'@md

# GET FITTED VALUES
# WRITE INITS FUNCTIONS FOR EACH PRIOR FAMILY

# occupancy = ~ psi_cov1
# detection = ~ p_cov1
# niter = 1000
# nchains = 3
# seed = NULL
# save_model = FALSE
# model_name = paste0("occ_model_", Sys.Date())

occ_mod <- function(occupancy, detection, data, niter = 1000, nchains = 3, seed = NULL,
                    save_model = FALSE, model_name = paste0("occ_model_", Sys.Date()),
                    beta_prior = "dnorm(0, 1/2)"){
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
    "for(i in 1:", ncol(occupancy_mod_matrix), "){beta_psi[i] ~ ", beta_prior, "}",
    "\n",
    "for(i in 1:", ncol(detection_mod_matrix), "){beta_p[i] ~ ", beta_prior, "}",
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
  message("\nDone.\n")

  # compute fitted values
  message("\nComputing fitted values.\n")
  for(i in 1:nchains){
    # psi
    betas <- samples[[i]][,grepl("_psi[[]", colnames(samples[[i]]))]
    tmp <- exp(betas %*% t(occupancy_mod_matrix))
    psi <- tmp / (1 + tmp)
  }
  cat("\nDone.\n")

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

#' Print method for `occ_mod` class
#'
#' @param x An object of class [occ_mod]
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

#' Summary method for `occ_mod` class
#'
#' @param x An object of class [occ_mod]
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
  n_eff <- coda::effectiveSize(coda::as.mcmc.list(lapply(reg_coefs, function(x) coda::as.mcmc(x))))

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
    `0.975` = apply(occupancy_coef, 2, quantile, prob = .975),
    n_eff = n_eff[grepl("_psi[[]", names(n_eff))]
  )

  # detection coefficients
  detection_coef <- samples_burnin[, grepl("_p[[]", colnames(samples_burnin))]
  detect_coef <- cbind(
    estimate = colMeans(detection_coef),
    sd = apply(detection_coef, 2, sd),
    `0.025` = apply(detection_coef, 2, quantile, prob = .025),
    `0.5` = apply(detection_coef, 2, quantile, prob = .5),
    `0.975` = apply(detection_coef, 2, quantile, prob = .975),
    n_eff = n_eff[grepl("_p[[]", names(n_eff))]
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
  attr(out, "nchains") <- nchains
  attr(out, "niter") <- niter
  attr(out, "burnin") <- burnin
  return(out)
}

#' Print method for `summary.occ_mod` class
#'
#' @param x An object of class [summary.occ_mod]
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

  # mcmc info
  cat(
    paste0(
      "\nMCMC information:\n",
      "nchains: ", attr(x, "nchains"), "\n",
      "niter: ", attr(x, "niter"), "\n",
      "burnin: ", attr(x, "burnin"), "\n"
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



