# simulate data
n.species <- 4
sim <- sim_MSMocc(K = n.species, M = 100, M.conf = 10, max_j = 20,
                 beta_psi = matrix(c(0, 1,
                                     0, 1.1,
                                     0, 1.2,
                                     0, 1.3),
                                   ncol = n.species),
                 beta_lambda = matrix(c(0, 1, 1,
                                        0, 1, 1,
                                        0, 1, 1.2,
                                        0, 1, 1.4),
                                      ncol = n.species),
                 theta = matrix(c(0.7, 0.1, 0.1, 0.1,
                                  0.1, 0.7, 0.1, 0.1,
                                  0.1, 0.1, 0.7, 0.1,
                                  0.1, 0.1, 0.1, 0.7), nrow = n.species),
                 rand_visits = TRUE)
data <- sim$data

# note structure of data frame
names(data)
head(data)

# fit model
ex <- MSMocc_mod(occupancy = ~ psi_cov1,
                 detection = ~ lambda_cov1 + lambda_cov2,
                 K = 4, species.names = names(data)[3:6],
                 data = data, niter = 40, beta_prior = "dnorm(0,2)")

# results
ex
summary(ex)
str(ex)
