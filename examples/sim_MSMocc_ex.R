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
sim$Y
