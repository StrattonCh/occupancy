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
