#'  Simulate data from the multi-species, misclassification count model for a single season.
#'
#' @description This function simulates data from he multi-species, misclassification count occupancy model first developed by
#'  \href{https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13315}{Wright et al. (2020)}.
#'
#'@details This function simulates data from the single season, multiple
#'  species occupancy model to estimate occupancy and relative activity for multiple species simultaneously. If `rand_visits
#'  = TRUE`, each site is visited a random number of times between two and
#'  `max_j`. Covariates are drawn from the uniform(0, 1.5) distribution so
#'  that the effect of the direction of each regression coefficient is
#'  intuitive. Note that if covariates are not desired, `beta_psi` and
#'  `beta_lambda` can be set to intercepts that generate the desired derived
#'  probabilities.
#'
#' @param K number of species
#' @param M number of sites
#' @param M.conf number of sites that have confirmed observations
#' @param max_j maximum number of visits to each site. If `rand_visits =
#'  FALSE`, this value is the number of visits to each site. If
#'  `rand_visits = TRUE`, each site is visited a random
#'  (`sample(2:max_j, size = 1)`) number of times.
#' @param beta_psi vector of regression coefficients used to generate psi
#' @param beta_lambda vector of regression coefficients used to generate lambda
#' @param seed optional seed for reproducibility
#'@param rand_visits logical; should each site be visited a random number of
#'  times? See details.
#' @param theta confusion matrix for automatic call classifications. Rows of theta should sum to one, diagonal elements are correct identification probabilities, first row corresponds to probabilities for calls truly from species 1.
#'
#' @return object of class `list` containing the following elements: \cr
#'  * `beta_psi` vector of regression coefficients used to generate psi
#'  * `beta_lambda` vector of regression coefficients used to generate lambda
#'  * `psi_cov` matrix of site level covariates
#'  * `lambda_cov` array of detection level covariates; each slice represents a
#'  single covariate
#'  * `psi` vector of derived site level occupancy probabilities
#'  * `lambda` matrix of derived visit level expected count detection parameter
#'  * `z` matrix of latent occupancy states for each site/species
#'  * `Y` array of observed Bernoulli responses for site/visit/species
#'  * `n_visits` vector of number of visits to each site
#'  * `c.conf` how calls confirmed to species K were classified for confirmed visits
#'  * `c.sum.a` and `c.sum.b`  unconfirmed calls classified TO species K from sites with some confirmations ("a") or no confirmations ("b")
#'  * `sites.conf` site locations with confirmations
#'  * `sites.unconf.a` and `sites.unconf.b` site locations without confirmations ("a" indecates some observations at these sites were confirmed, "b" sites had no confirmations)
#'  * `data` a data frame containing all information necessary to fit the
#'  model
#' @export
#'
#' @importFrom stats runif rbinom rmultinom
#'
#' @example examples/sim_MSMocc_ex.R
#'
sim_MSMocc <- function(K = 2, M = 20, M.conf = 2, max_j = 10,
                    beta_psi = matrix(c(0, 1, 0, 1.2), ncol = K),
                    beta_lambda = matrix(c(0, 1, 0, 1.2), ncol = K),
                    theta = matrix(c(0.8, 0.2, 0.2, 0.8), byrow = T, ncol = K),
                    seed = NULL, rand_visits = TRUE) {

  # check that number sites with confirmed observations is less than total number of sites
  if (M.conf > M) stop("The number of sites with confirmed observations (M.conf) must be less than or equal to the total number of sites (M).")

  # create out vector
  out <- list()

  # optional seed
  if (!is.null(seed)) set.seed(seed)

  # convenience
  inv.logit <- function(x) exp(x) / (1 + exp(x))

  # double check betas; separate into list by species
    beta_psi <- as.matrix(beta_psi, ncol = K)
    beta_lambda <- as.matrix(beta_lambda, ncol = K)

  # dimension of regression coefficients
  # assumption: all species have same number of coefficients
    p_beta_psi <- nrow(beta_psi)
    p_beta_lambda <- nrow(beta_lambda)

  # generate covariates
  ## psi
    if (p_beta_psi == 1) {
      psi_cov <- matrix(rep(1, M), ncol = 1)
    } else{
      psi_cov <- cbind(
        rep(1, M),
        matrix(
          stats::runif((p_beta_psi - 1) * M, 0, 1.5), ncol = p_beta_psi - 1
        )
      )
    }

  # column names
  if (p_beta_psi == 1) {
    colnames(psi_cov) <- "psi_int"
  } else {
    colnames(psi_cov) <- c("psi_int", paste0("psi_cov", 1:(ncol(psi_cov) - 1)))
  }

  # generate psi for each species (one column per species)
    psi <- inv.logit(psi_cov %*% beta_psi)

  # lambda
  lambda_cov <- array(0, dim = c(M, max_j, p_beta_lambda))
  lambda_linpred <- array(0, dim = c(M, max_j, K))
  for (k in 1:p_beta_lambda) {
    if (k == 1) {
      lambda_cov[, , k] <- 1
    } else {
      lambda_cov[, , k] <- stats::runif(n = M * max_j, 0, 1.5)
    }
    for (i in 1:K) {
      lambda_linpred[, , i] <- lambda_linpred[,, i] + beta_lambda[k, i] * lambda_cov[, , k]
    }
  }

  # names
  if (p_beta_lambda == 1) {
    dimnames(lambda_cov)[[3]] <- "lambda_int"
  } else {
    dimnames(lambda_cov)[[3]] <- c("lambda_int", paste0("lambda_cov", 1:(dim(lambda_cov)[3] - 1)))
  }

  # generate lambda
  lambda <- array(0, dim = c(M, max_j, K))
  for (i in 1:K) {
    lambda[ , , i] <- inv.logit(lambda_linpred[, , i])
  }

  # generate latent occupancy for each species
  # z <- list()
  # for (i in 1:K) {
    z <- apply(psi, c(1,2), function(x) stats::rbinom(1, 1, x))
  # }

  # generate responses
  Y <- array(0, dim = c(M, max_j, K))
  for (k in 1:K) {
    for (i in 1:M) {
      for (j in 1:max_j) {
        Y[i, j, k] <- z[i,k] * stats::rpois(1, lambda = lambda[i, j, k])
      }
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
        Y[i, (n_visits[i] + 1):max_j, ] <- NA

        # covariates
        lambda_cov[i, (n_visits[i] + 1):max_j, ] <- NA

        # derived parameters
        lambda[i, (n_visits[i] + 1):max_j, ] <- NA
      }
    }

    out$n_visits <- apply(Y, 1, function(x) length(which(!is.na(x)))/K)
  }
  if (!rand_visits) out$n_visits <- rep(max_j, M)

  # generate classifications

  c <- array(NA, dim = c(M * max_j, K, K))
  for (k in 1:K) {
      c[, ,k] <- t(sapply(as.vector(Y[, , k]),
                          function(x) {
                            if (is.na(x) == T) {
                              return(rep(NA, K))
                            }else {rmultinom(1, x, theta[k, ])}
                            }))
  }

  #total number of observations
  n.obs <-  sum(out$n_visits)
  #total number of observations at sites with (some) confirmed observations
  n.obs.conf <- sum(out$n_visits[1:M.conf])
  #number of confirmed observations per site with confirmed observations
  n.conf <- 2
  n.unconf <- out$n_visits[1:M.conf] - n.conf
  keep.conf <- vector()

    for (i in 1:M.conf) {
     keep.conf <- c(keep.conf,  c(rep(TRUE, n.conf), rep(FALSE, n.unconf[i])))
    }

  # c.conf[obs, species (how classified), species (which confirmed)] is how calls confirmed to each species were classified for confirmed visits
  c.conf <- c[1:n.obs.conf, , ][keep.conf, ,]
  # same as above but for unconfirmed visits
  c.unconf <- c[1:n.obs.conf, ,][!keep.conf, ,]

  # c.sum.a[unconfirmed calls (at sites with confirmed calls), species] is the unconfirmed calls classified TO species k from sites that had some confirmations
  # c.sum.b is counts of unconfirmed calls to species k from sites with no confirmations
  c.sum.a <- matrix(NA, ncol = K, nrow = dim(c.unconf)[1])
  c.sum.b <- matrix(NA, ncol = K, nrow = length((n.obs.conf + 1):n.obs))
  for (i in 1:K) {
    c.sum.a[, i] <- apply(c.unconf[, i, ], 1, sum, na.rm = T)

    c.sum.b[, i] <-  apply(c[(n.obs.conf + 1):n.obs, i, ], 1, sum, na.rm = T)
  }

  # all site variables just keep track of which sites correspond to each set of observations

  sites.conf <- vector()
  for (i in 1:M.conf) {
    sites.conf <- sites.unconf <- c(sites.conf, rep(i, out$n_visits[i]))
  }
  sites.conf <- sites.conf[keep.conf]
  sites.unconf.a <- sites.unconf[!keep.conf]
  sites.unconf.b <- vector()
  for (i in (M.conf + 1):M) {
    sites.unconf.b <- c(sites.unconf.b, rep(i, out$n_visits[i]))
  }

  # return
  out$beta_psi <- beta_psi
  out$beta_lambda <- beta_lambda
  out$psi_cov <- psi_cov
  out$lambda_cov <- lambda_cov
  out$psi <- psi
  out$lambda <- lambda
  out$z <- z
  out$Y <- Y
  out$c.conf = c.conf
  out$c.sum.a = c.sum.a
  out$c.sum.b = c.sum.b
  out$sites.conf = sites.conf
  out$sites.unconf.a = sites.unconf.a
  out$sites.unconf.b = sites.unconf.b

  # create data compatible with occ_mod

  # data <- list(
    data = data.frame(
      site = rep(1:M, out$n_visits),
      visit = unlist(c(sapply(out$n_visits, function(x) 1:x)))
    )
  # )

  # multiple species data
  for (i in 1:K) {
    data <- cbind(data, stats::na.omit(c(t(out$Y[,,i]))))
  }
  names(data)[-c(1:2)] <- paste("species", 1:K)

  # site covariates
  tmp <- out$psi_cov
  tmp[,1] <- 1:nrow(tmp);colnames(tmp)[1] <- "site"
  tmp <- as.data.frame(tmp)
  data <- dplyr::left_join(data, tmp, by = "site")

  # detection covariates
  tmp <- as.data.frame(apply(out$lambda_cov, 3, function(x) stats::na.omit(c(t(x)))))
  tmp$site <- data$site
  tmp$visit <- data$visit
  tmp <- tmp[,-which(names(tmp) == "lambda_int")]
  data <- dplyr::left_join(data, tmp, by = c("site", "visit"))

  # export
  out$data <- data

  return(out)
}


#' Title
#'
#' @param occupancy
#' @param detection
#' @param data
#' @param nspecies
#' @param niter
#' @param nchains
#' @param seed
#' @param save_model
#' @param model_name
#'
#' @return
#' @export
#'
#' @examples
#'
#' @importFrom rlang .data
#'
MSMocc_mod <- function(occupancy, detection, data, K, species.names,
                    niter = 1000, nchains = 3, seed = NULL,
                    save_model = FALSE, model_name = paste0("MSMocc_model_", Sys.Date()),
                    beta_prior = "dnorm(0, 1/2)"){

  # convenience
  `%notin%` <- Negate("%in%")

  # check names on data
  if ("site" %notin% names(data)) stop("Data must contain column named 'site'")
  if ("visit" %notin% names(data)) stop("Data must contain column named 'visit'")
  # if ("y" %notin% names(data)) stop("Data must contain column named 'y'")
  for (i in 1:K) {
    if (species.names[i] %notin% names(data)) stop("Data must contain a column for each species matching provided species names.")
  }


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
  if (nrow(occupancy_mod_matrix) != length(unique(data$site))) stop("There are more unique site-level covariate values than there are sites.")
  if (occ_num_params == 1) {
    occupancy_mod_matrix <- as.matrix(occupancy_mod_matrix[,-which(colnames(occupancy_mod_matrix) == "site_")], ncol = 1)
    colnames(occupancy_mod_matrix) <- "(Intercept)"
  } else{
    occupancy_mod_matrix <- occupancy_mod_matrix[,-which(colnames(occupancy_mod_matrix) == "site_")]
  }

  ## detection
  detection_mod_matrix <- stats::model.matrix(object = detection_mod, data = data)

  ## number of sites (M) and visits (n_visits)
  M <- length(unique(data$site))
  n_visits <- data %>%
    dplyr::select(.data$site, .data$visit) %>%
    dplyr::group_by(.data$site) %>%
    dplyr::arrange(dplyr::desc(.data$visit)) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::pull(.data$visit) %>%
    as.vector()

  # prep data for NIMBLE
  nimble_data <- list()

  # prep list for multiple species' data
  Y <- array(NA, c(M, max(n_visits), K))
  ## responses
  for (i in 1:K) {
    Y[,,i] <- data %>%
      dplyr::select(.data$site, .data$visit, .data[[species.names[i]]]) %>%
      dplyr::group_by(.data$site) %>%
      tidyr::pivot_wider(names_from = "visit",
                         values_from = species.names[i]) %>%
      dplyr::ungroup() %>%
      dplyr::select(`1`:ncol(.)) %>%
      as.matrix()
  }

  nimble_data$Y <- Y

  ## occupancy level
  occupancy_covariates <- list()
  for (i in 1:ncol(occupancy_mod_matrix)) {
    occupancy_covariates[[colnames(occupancy_mod_matrix)[i]]] <- unname(occupancy_mod_matrix[,i])
  }
  nimble_data$occupancy_covariates <- occupancy_covariates

  ## detection level
  detection_covariates <- list()
  for (i in 1:ncol(detection_mod_matrix)) {
    detection_covariates[[colnames(detection_mod_matrix)[i]]] <- data %>%
      dplyr::select(.data$site, .data$visit) %>%
      dplyr::mutate(tmp = detection_mod_matrix[,i]) %>%
      dplyr::group_by(.data$site) %>%
      tidyr::pivot_wider(names_from = "visit", values_from = "tmp") %>%
      dplyr::ungroup() %>%
      dplyr::select(`1`:ncol(.)) %>%
      as.matrix()
  }
  nimble_data$detection_covariates <- detection_covariates

  ## nvisits

   nimble_data$nvisits <- n_visits

   ###### reorganize sites and visits into observations to loop through ########
   obs <- 1:(sum(n_visits))
   nimble_data$obs <- obs

   # NIMBLE code
   ## priors
   priors <- paste0(
     "# priors\n",
     "for(i in 1:", K, "){ \n",
     "for(j in 1:", ncol(occupancy_mod_matrix), "){beta_psi[i,j] ~ ", beta_prior, "}",
     "\n",
     "for(j in 1:", ncol(detection_mod_matrix), "){beta_lambda[i,j] ~ ", beta_prior, "}",
     "}",
     "\n"
     ## should theta be a provided matrix or have a prior??
   )

   #--------------------------------------------------------------
   # likelihoods
   #--------------------------------------------------------------
   likelihood_a <- paste0(
     "# likelihood \n",
     "# loop through species and sites \n",
     "for(i in 1:",K,"){ \n",
     "for(j in 1:",M,"){ \n",
     " # calculate psi \n")
   if (ncol(occupancy_mod_matrix) > 1) {
     likelihood_b <- paste0(
       " logit(psi[i,j]) <- beta_psi[i, 1] +\n",
       paste0(sapply(2:ncol(occupancy_mod_matrix), function(x){
         paste0("  beta_psi[i,", x, "] * ", colnames(occupancy_mod_matrix)[x], "[j]")
       }), collapse = " + \n")
     )
   } else{
     likelihood_b <- paste0(
       " logit(psi[i,j]) <- beta_psi[i,1]"
     )
   }

   likelihood_c <- paste0(
     "\n\n",
     " # latent occupancy state \n",
     " z[i,j] ~ dbinom(psi[i,j], 1) \n\n",
     " # loop through visits \n",
     " for(k in 1:nvisits[j]){\n",
     "  # calculate lambda\n"
   )
   if (ncol(detection_mod_matrix) > 1) {
     likelihood_d <- paste0(
       "  logit(lambda[i,j,k]) <- beta_lambda[i,1] +\n",
       paste0(sapply(2:ncol(detection_mod_matrix), function(x){
         paste0("   beta_lambda[i,", x, "] * ", colnames(detection_mod_matrix)[x], "[j,k]")
       }), collapse = " + \n")
     )
   } else{
     likelihood_d <- paste0(
       "  logit(p[i,j,k]) <- beta_lambda[i,1]"
     )
   }
   likelihood_e <- paste0(
     "\n\n",
     "  # response model\n",
     "  y.true[j,k, i] ~ dpois(lambda[i,j,k]*z[i,j])",
     "\n }",
     "\n} \n",
     "y.true.obs[1:",obs, ", i] <- y.true[1:", M, ", 1:", n_visits, ", i] \n",
     "y.counts[i] <- sum(y.true[1:M,1:max_visits, i])",
     "\n}"
   )

   likelihood_f <- paste0(
     "\n\n",
     " #classification observation model\n",
     "for (k in 1:",K, ") {\n",
       "for (obs in 1:",obs, ") { \n",

           "C.obs[obs,1:",K," ,k] ~ dmulti(theta[k, 1:",K, "], y.true.obs[obs, k])\n",
     "\n}",
     "\n}"
   )

   likelihood <- paste0(likelihood_a, likelihood_b, likelihood_c, likelihood_d, likelihood_e, likelihood_f)

   code <- paste0(priors, "\n", likelihood)

   # initialization
   mod_inits <- function(x){
     tmp <- x
     list(
       z = t(apply(nimble_data$Y[, ,], c(1, 3),
                   function(x) ifelse(sum(x, na.rm = T) == 0, 0, 1))),
       beta_psi = matrix(stats::rnorm(K * ncol(occupancy_mod_matrix), 0, sqrt(2)),
                         nrow = K),
       beta_lambda = matrix(stats::rnorm(K * ncol(detection_mod_matrix), 0, sqrt(2)),
                            nrow = K),
       y.counts = apply(nimble_data$Y[, , ], c(3),
                        sum, na.rm = T)
     )
   }

   # model
   data_ <- c(nimble_data$occupancy_covariates, nimble_data$detection_covariates)
   data_ <- data_[-c(which(names(data_) == "(Intercept)"))]
   data_[["M"]] <- nrow(nimble_data$Y)
   data_[["K"]] <- dim(nimble_data$Y)[3]
   data_[["nvisits"]] <- nimble_data$nvisits
   data_[["max_visits"]] <- max(nimble_data$nvisits)
   data_[["Y"]] <- nimble_data$Y

   # write model
   if (save_model) message("Writing NIMBLE model to current working directory:")
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
   if (!save_model) file.remove(paste0(model_name, ".txt"))

   # fit model
   message("\nCompiling R model to C")
   model_c <- nimble::compileNimble(model)

   message("\nBuilding R MCMC object")
   mcmcinfo <- utils::capture.output(mcmc <- nimble::buildMCMC(model_c,
                                                               monitors = c("z",
                                                                            "beta_psi",
                                                                            "beta_lambda")))

   message("\nCompiling R MCMC object to C")
   mcmc_c <- nimble::compileNimble(mcmc)

   message("\nRunning MCMC")
   inits <- lapply(1:nchains, function(x) mod_inits(x))
   samples <- nimble::runMCMC(mcmc_c, niter = niter, inits = inits, nchains = nchains)

   # compute posterior likelihood
   # compute posterior likelihood
   message("\nComputing posterior likelihood")

   multispecies_loglik <- list()

   for (k in 1:K) {
     multispecies_loglik[[species.names[k]]] <- list()
     for (i in 1:nchains) {
       samples_ <- samples[[i]]
       beta_psi <- samples_[,which(colnames(samples_) %in% paste0('beta_psi[',k, ", ", 1:ncol(occupancy_mod_matrix), ']'))]
       beta_lambda <- samples_[,which(colnames(samples_) %in% paste0('beta_lambda[',k, ", ", 1:ncol(detection_mod_matrix), ']'))]
       z <- samples_[,which(colnames(samples_) %in% paste0('z[',k, ", ", 1:nrow(data_[["Y"]]), ']'))]

       z_ <- z[,rep(1:ncol(z),data_[["nvisits"]])]
       y_ <- matrix(rep(data[, species.names[k]],  nrow(samples_)), nrow = nrow(samples_), byrow = TRUE)

       psi_logit <- exp(beta_psi %*% t(occupancy_mod_matrix))
       psi <- psi_logit / (1 + psi_logit)
       lambda_logit <- exp(beta_lambda %*% t(detection_mod_matrix))
       lambda <- lambda_logit / (1 + lambda_logit)

       # log pointwise predictive density
       site <- psi^z * (1 - psi)^(1 - z)

       visit <- ((z_*lambda)^y_ * exp^(-z_ * lambda))/(factorial(y_))
       ndx_upr <- cumsum(data_[["nvisits"]])
       ndx_lwr <- ndx_upr - data_[["nvisits"]] + 1
       tmp <- sapply(1:ncol(site), function(x){ apply(visit[,ndx_lwr[x]:ndx_upr[x]], 1, prod)})

       multispecies_loglik[[species.names[k]]][[i]] <- log(site * tmp)
       colnames( multispecies_loglik[[species.names[k]]][[i]]) <- paste0("site", 1:nrow(data_[[species.names[k]]]))
     }
   }
   names( multispecies_loglik[[species.names[k]]]) <- paste0("chain", 1:nchains)
   message("\nDone.\n")

   # compute fitted values
   message("\nComputing fitted values.\n")
   psi_mcmc <- list()
   lambda_mcmc <- list()
   for (i in 1:nchains) {
     # psi
     betas <- samples[[i]][,grepl("_psi[[]", colnames(samples[[i]]))]
     tmp <- exp(betas %*% t(occupancy_mod_matrix))
     psi <- tmp / (1 + tmp)
     colnames(psi) <- paste0("site", 1:ncol(psi))
     psi_mcmc[[paste0("chain", i)]] <- psi

     # detection
     betas <- samples[[i]][,grepl("_lambda[[]", colnames(samples[[i]]))]
     lambda_mcmc_chain <- array(0, dim = c(
       nrow(psi), nrow(detection_covariates[[1]]), ncol(detection_covariates[[2]])
     ))
     for (iter in 1:nrow(psi)) {
       beta_iter <- betas[iter,]
       tmp <- sapply(1:length(detection_covariates), function(x) beta_iter[x] * detection_covariates[[x]], simplify = "array")
       tmp2 <- exp(apply(tmp, 2, rowSums))
       lambda_mcmc_chain[iter,,] <- tmp2 / (1 + tmp2)
       dimnames(lambda_mcmc_chain)[[2]] <- paste0("site", 1:(dim(lambda_mcmc_chain)[2]))
       dimnames(lambda_mcmc_chain)[[3]] <- paste0("visit", 1:(dim(lambda_mcmc_chain)[3]))

     }
     lambda_mcmc[[paste0("chain", i)]] <- lambda_mcmc_chain

   }
   cat("\nDone.\n")

   # return
   out <- list(samples = samples, loglik = loglik,
               derived = list(psi = psi_mcmc,
                              lambda = lambda_mcmc))

   # attributes
   class(out) <- c("MCMocc_mod", "list")
   attr(out, "code") <- code
   attr(out, "mcmcinfo") <- mcmcinfo
   attr(out, "occupancy_call") <- occupancy
   attr(out, "detection_call") <- detection
   return(out)

     }


