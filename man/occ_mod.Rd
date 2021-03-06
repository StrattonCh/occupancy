% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/occ_functions.R
\name{occ_mod}
\alias{occ_mod}
\title{Fit the single species, single season site occupancy model.}
\usage{
occ_mod(
  occupancy,
  detection,
  data,
  niter = 1000,
  nchains = 3,
  seed = NULL,
  save_model = FALSE,
  model_name = paste0("occ_model_", Sys.Date()),
  beta_prior = "dnorm(0, 1/2)"
)
}
\arguments{
\item{occupancy}{model declaration for the occupancy portion of the model
using standard linear model syntax}

\item{detection}{model declaration for the detection portion of the model
using standard linear model syntax}

\item{data}{data from which to create the detection and covariate matrices.
Each row should represent a single visit to a site. See details and examples
for proper formatting}

\item{niter}{number of MCMC iterations}

\item{nchains}{number of MCMC chains}

\item{seed}{optional seed}

\item{save_model}{logical; should a text file containing the model be
exported?}

\item{model_name}{character string defining the name of the text file
describing the model if \code{save_model = TRUE}}

\item{beta_prior}{character string defining prior distribution for regression
coefficients at the occupancy and detection levels. Priors should be
specified using distributions available in NIMBLE. See
\href{https://r-nimble.org/html_manual/cha-writing-models.html#subsec:dists-and-functions). et al. (2002)}{available distributons in NIMBLE}}
}
\value{
an object of class \code{list} containing the following:\cr
\itemize{
\item \code{samples} object of class \code{list} of length \code{nchains}, each
containing a \code{matrix} posterior samples
\item \code{loglik} object of class \code{list} of length code{nchains}, each
containing a \code{matrix} of samples of the log posterior likelihood
}
}
\description{
This function fits the single species, single
season site occupancy model first developed by
\href{https://esajournals.onlinelibrary.wiley.com/doi/10.1890/0012-9658\%282002\%29083\%5B2248\%3AESORWD\%5D2.0.CO\%3B2}{MacKenzie et al. (2002)}.
}
\details{
This function fits the single season, single species site occupancy
model using the logit link function. The \code{data} should contain columns
named \code{site}, \code{visit}, and \code{y}. See examples.
}
\examples{
# simulate data
sim <- sim_occ(M = 100, max_j = 10, seed = 01012021, rand_visits = FALSE)
data <- sim$data

# note structure of data frame
names(data)
head(data)

# fit model
ex <- occ_mod(occupancy = ~ psi_cov1, detection = ~ p_cov1, data = data,
              niter = 4000, beta_prior = "dunif(-5, -2)")

# results
ex
summary(ex)
str(ex)
}
