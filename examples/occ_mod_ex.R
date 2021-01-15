# simulate data - note structure of data
sim <- sim_occ(M = 250, max_j = 10, seed = 11201995, rand_visits = FALSE)
data <- sim$data
names(data)

# fit model
ex <- occ_mod(occupancy = ~ psi_cov1, detection = ~ p_cov1, data = data,
              niter = 5000)
colMeans(ex)
plot(ex[,4], type = "l")
