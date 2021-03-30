library(mice) # data
library(INLA)
library(mvtnorm)

# loading general functions
source("./genFuncs.R")
# loading amis with inla functions
source("./inlaAMIS.R")
# loading is with inla functions
source("./inlaIS.R")
# loading mcmc with inla functions
source("./inlaMH.R")


data(nhanes2)

d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)

df = list(d.mis = d.mis, idx.mis = idx.mis)

init = list(mu = rep(mean(df$d.mis$bmi, na.rm = TRUE),n.mis),
            cov = diag(2*var(df$d.mis$bmi, na.rm = TRUE),n.mis,n.mis))

names(init$mu) = sprintf("Observation_%d",df$idx.mis)
colnames(init$cov) = sprintf("Observation_%d",df$idx.mis)

fit.inla <- function(data, eta) {

  data$d.mis$bmi[data$idx.mis] = eta

  res = inla(chl ~ 1 + bmi + age, data = data$d.mis)

  return(list(mlik = res$mlik[[1]],
              dists = list(beta0 = res$marginals.fixed[[1]],
                           beta1 = res$marginals.fixed[[2]],
                           beta2 = res$marginals.fixed[[3]],
                           beta3 = res$marginals.fixed[[4]],
                           tau = res$marginals.hyperpar[[1]])))
}

prior.x.mis <- function(x, mu = mean(d.mis$bmi, na.rm = TRUE),
                        sigma = 2*sd(d.mis$bmi, na.rm = TRUE), log = TRUE) {
  res = dnorm(x, mean = mu, sd= sigma, log = log)
  if(log) {
    return(sum(res))
  } else {
    return(prod(res))
  }
}

mc
dq.x.mis <- function(y, x, sigma = sqrt(5), log =TRUE) {
  sum(dnorm(x, mean = y, sd = sqrt(5), log = log))
}

rq.x.mis <- function(x, sigma = sqrt(5)) {
  rnorm(length(x), mean = x, sd = sqrt(5))
}

### AMIS
amis_mod <- inlaAMIS(data = df, init = init, prior.x.mis,
                               dq.x.mis, rq.x.mis, fit.inla,
                              N_t = seq(25,50,1)*10, N_0 = 250,ncores= 10)
save(amis_mod, file = "./missing//missing-amis.Rdata")

### IS
is_mod <- inlaIS(data = df, init = init, prior.x.mis,
                           dq.x.mis, rq.x.mis,fit.inla, N_0 = 800, N = 10000,ncores = 10)
save(is_mod, file = "./missing//missing-is.Rdata")

### MCMC
mcmc_mod <- inlaMH(data = df, init =init,
                               prior.x.mis, dq.x.mis, rq.x.mis, fit.inla,
                               n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_mod, file = "./missing//missing-mcmc.Rdata")


