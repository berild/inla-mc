# Fit ZIP models with MCMC using jags
# --"-- with INLA-IS
# --"-- with INLA-AMIS

library(INLA)
library(pscl)
library(mvtnorm)
# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("./inlaMC/inlaMC.R")

# Load data
# https://stats.idre.ucla.edu/r/dae/zip/
#zinb <- read.csv("https://stats.idre.ucla.edu/stat/data/fish.csv")
zinb <- read.csv("./zip/fish.csv")
zinb <- within(zinb, {
    nofish <- factor(nofish)
    livebait <- factor(livebait)
    camper <- factor(camper)
})

summary(zinb)


# Model:
#  Poisson: count ~ child + camper 
# Binomial: ~ 1 + persons

#summary(m1 <- zeroinfl(count ~ child + camper | persons, data = zinb))


# Reduce data for testing
#d <- zinb[1:30, ]

# Full data
d <- zinb

# ZIP
ml.res <- summary(zeroinfl(count ~ child + camper | persons, data = d))
ml.res

# Fit model with JAGS
library(rjags)
d.jags <- list(N = nrow(d), count = d$count, persons = d$persons,
  child = d$child, camper = d$camper)

# Set seed
inits.jags <- list(".RNG.name"="base::Mersenne-Twister",
  ".RNG.seed" = 1)

jm <- jags.model (file = "./zip/ZIP.bug", data = d.jags, inits = inits.jags)

# Burn-in
update(jm, 10000)

js <- jags.samples(jm, n.iter = 50000, thin = 10,
  variable.names = c("gamma", "beta", "p", "mu"))

# Summary
js$gamma
js$beta

# Gamma
apply(js$gamma[,,1], 1, mean)
apply(js$gamma[,,1], 1, sd)

# Beta
apply(js$beta[,,1], 1, mean)
apply(js$beta[,,1], 1, sd)

#save
save(file = "./sims/zip/mcmc_zip.RData", list = c("jm", "js", "d.jags", "inits.jags"))

# Fit model with IS/AMIS with INLA
fit.inla <- function(data, eta) {

  logit_pi <- eta[1] + eta[2] * data$persons

  # Define pi hyper for likelihood
  hyper_pi <- lapply(logit_pi, function(X) {
    list(hyper = list(prob = list(fixed = TRUE, initial = X)))
  })

  # Define COUNT as diagonal matrix of observed counts
  COUNT <- matrix(NA, nrow = nrow(data), ncol = nrow(data))
  diag(COUNT) <- data$count

  res <- inla(COUNT ~ child + camper, 
    data = list(COUNT = COUNT, child = d$child, camper = d$camper),
    family = rep("zeroinflatedpoisson1", nrow(d)),
    num.threads = "1:1",
    control.fixed = list(prec.intercept = 0.001),
    control.family = hyper_pi
  )
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           child = res$marginals.fixed[[2]],
                           camper1 = res$marginals.fixed[[3]])))
}

# Prior for 
prior.beta <- function(x, y = 0, sigma = sqrt(1000), log = TRUE) {
  sum(dnorm(x = x, mean = y, sd = sigma, log = TRUE))
}

## Proposal distribution

# Moments of proposal
# NOTE: We use ML estimates to set the sampling distribution
# This is fine, in principle, becuase then samples are re-weighted
init = list(mu = as.numeric(ml.res$coefficients[[2]][1:2,1]),cov = diag(3*as.numeric(ml.res$coefficients[[2]][1:2,2])))
init

# Evaluate
dq.beta <- function(y, theta = init, log =TRUE) {
  dmvt(y,delta=theta[[1]],sigma=theta[[2]],df=3,type = "shifted",log = log) 
}
# Sample
rq.beta <- function(theta = init) {
  as.vector(rmvt(1,sigma = theta[[2]], df=3, delta = theta[[1]], type = "shifted"))
}


## Fit models
n.samples = 10000
# IS-INLA
set.seed(1)
is_mod <- inlaIS(data=d, init=init, prior = prior.beta, d.prop=dq.beta, r.prop=rq.beta, fit.inla=fit.inla, N = n.samples, kde = TRUE, ncores=60)
save(is_mod, file = "./sims/zip/is_zip.Rdata")
# AMIS-INLA
set.seed(1)
amis_mod <- inlaAMIS(data=d, init=init, prior= prior.beta, d.prop=dq.beta, r.prop=rq.beta, fit.inla=fit.inla, N_t = rep(2000, 4), N_0 = 2000, kde = TRUE, ncores = 60)
save(amis_mod, file = "./sims/zip/amis_zip.Rdata")
