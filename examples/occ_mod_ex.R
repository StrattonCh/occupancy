# simulate data
sim <- sim_occ(M = 250, max_j = 10, seed = 01012021, rand_visits = FALSE)
data <- sim$data

# note structure of data frame
names(data)
head(data)

# fit model
ex <- occ_mod(occupancy = ~ psi_cov1, detection = ~ p_cov1, data = data,
              niter = 2000)

# results
str(ex)
ex



# lppd <- sum(log(colMeans(exp(loglik))))
# pwaic1 <- 2 * sum(log(colMeans(exp(loglik))) - colMeans(loglik))
# pwaic2 <- sum(apply(loglik, 2, var))
# waic1 <- -2 * (lppd - pwaic1)
# waic2 <- -2 * (lppd - pwaic2)
