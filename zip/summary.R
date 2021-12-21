# Compute summary statistics in the zero-inflated poisson example
library(INLA)
library(rjags)

# JAGS
load("./sims/zip/mcmc_zip.RData")

# Summary
js$gamma
js$beta

# Gamma
apply(js$gamma[,,1], 1, mean)
apply(js$gamma[,,1], 1, sd)

# Beta
apply(js$beta[,,1], 1, mean)
apply(js$beta[,,1], 1, sd)


# IS
load("./sims/zip/is_zip.RData")

# Weigths
ww <-  is_mod$weight / sum( is_mod$weight)
summary(ww)

# n_e
1 / sum(ww^2)

# Gamma
lapply(is_mod$eta_kern, function(X) {unlist(inla.zmarginal(X, TRUE))})

# Beta
lapply(is_mod$margs, function(X) {unlist(inla.zmarginal(X, TRUE))})

# AMIS
load("./sims/zip/amis_zip.RData")

# Weigths
ww <-  amis_mod$weight / sum(amis_mod$weight)
summary(ww)

# n_e
1 / sum(ww^2)

# Gamma
lapply(amis_mod$eta_kern, function(X) {unlist(inla.zmarginal(X, TRUE))})

# Beta
lapply(amis_mod$margs, function(X) {unlist(inla.zmarginal(X, TRUE))})
