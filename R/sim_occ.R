
#'@title Simulate data from the single species, single season site occupancy model.
#'
#' @description This function simulates data from the single species, single
#'  season occupancy model first developed by
#'  \href{https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658%282002%29083%5B2248%3AESORWD%5D2.0.CO%3B2}{MacKenzie
#'   et al. (2002)}.
#'
#' @details This function simulates data from the vanilla single season, single
#'  species occupancy model using the logit link function. If \code{rand_visits
#'  = TRUE}, each site is visited a random number of times between two and
#'  \code{max_j}. Covariates are drawn from the uniform(0, 1.5) distribution so
#'  that the effect of the direction of each regression coefficient is
#'  intuitive. Note that if covariates are not desired, \code{beta_psi} and
#'  \code{beta_p} can be set to intercepts that generate the desired derived
#'  probabilities.
#'
#' @param M number of sites
#' @param max_j maximum number of visits to each site. If \code{rand_visits =
#'  FALSE}, this value is the number of visits to each site. If
#'  \code{rand_visits = TRUE}, each site is visited a random
#'  (\code{sample(2:max_j, size = 1)}) number of times.
#' @param beta_psi vector of regression coefficients used to generate psi
#' @param beta_p vector of regression coefficients used to generate p
#' @param seed optional seed for reproducibility
#' @param rand_visits logical; should each site be visited a random number of
#'  times? See details.
#'
#' @example examples/sim_occ_ex.R
#'
#' @return object of class \code{list} containing the following elements: \cr
#'  * \code{beta_psi} vector of regression coefficients used to generate psi
#'  * \code{beta_p} vector of regression coefficients used to generate p
#'  * \code{psi_cov} matrix of site level covariates
#'  * \code{p_cov} array of detection level covariates; each slice represents a
#'  single covariate
#'  * \code{psi} vector of derived site level occupancy probabilities
#'  * \code{p} matrix of derived visit level detection probabilities
#'  * \code{z} vector of latent occupancy states for each site
#'  * \code{Y} matrix of observed Bernoulli responses
#'  * \code{n_visits} vector of number of visits to each site
#'  * \code{data} a data frame containing all information necessary to fit the
#'  model
#'
#' @importFrom magrittr %>%
#' @export
#'
#' @md

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
  psi_cov <- cbind(
    rep(1, M),
    matrix(
      stats::runif((p_beta_psi - 1) * M, 0, 1.5),
      ncol = p_beta_psi - 1
    )
  )

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
    y = na.omit(c(t(out$Y)))
  )

  # site covariates
  tmp <- out$psi_cov
  tmp[,1] <- 1:nrow(tmp);colnames(tmp)[1] <- "site"
  tmp <- as.data.frame(tmp)
  data <- dplyr::left_join(data, tmp, by = "site")

  # detection covariates
  tmp <- as.data.frame(apply(out$p_cov, 3, function(x) na.omit(c(t(x)))))
  tmp$site <- data$site
  tmp$visit <- data$visit
  tmp <- tmp[,-which(names(tmp) == "p_int")]
  data <- dplyr::left_join(data, tmp, by = c("site", "visit"))

  # export
  out$data <- data

  return(out)
}
