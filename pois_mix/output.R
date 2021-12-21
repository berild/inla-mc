# Compute summary results
library("INLA")

# IS results
load("pois-mix/is_pois_mix.Rdata")

# Weights
ww <- is_mod$weight / sum(is_mod$weight)
summary(ww)

# n_e
1 / sum(ww^2)


# Summary results
# mu_1
inla.zmarginal(inla.tmarginal(exp, is_mod$margs[[1]]))
# mu_2
inla.zmarginal(inla.tmarginal(exp, is_mod$margs[[2]]))

# AMIS results
load("pois-mix/amis_pois_mix.Rdata")

# Weights
ww <- amis_mod$weight / sum(amis_mod$weight)
summary(ww)

# n_e
1 / sum(ww^2)


# Summary results
# mu_1
inla.zmarginal(inla.tmarginal(exp, amis_mod$margs[[1]]))
# mu_2
inla.zmarginal(inla.tmarginal(exp, amis_mod$margs[[2]]))


