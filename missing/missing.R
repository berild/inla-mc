library(mice) # data
library(INLA)
library(mvtnorm)

# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("./inlaMC/inlaMC.R")

# loading dataset
data(nhanes2)

# data
d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi)) # finding na's
n.mis <- length(idx.mis) # number of nans

df = list(d.mis = d.mis, idx.mis = idx.mis)

# initial parameters of the proposal distribution
init = list(mu = rep(mean(df$d.mis$bmi, na.rm = TRUE),n.mis),
            cov = diag(2*var(df$d.mis$bmi, na.rm = TRUE),n.mis,n.mis))
# Setting names on the imputed values
names(init$mu) = sprintf("Observation_%d",df$idx.mis)
colnames(init$cov) = sprintf("Observation_%d",df$idx.mis)

# conditional INLA call within the algorithms
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

# proposal distribution
## evaluate
dq.x.mis <- function(y, theta=init, log =TRUE) {
  dmvnorm(y, mean = theta[[1]], sigma = theta[[2]], log = log)
}
## sample
rq.x.mis <- function(theta=init) {
  as.vector(rmvnorm(1, mean = theta[[1]], sigma = theta[[2]]))
}

### AMIS
set.seed(1)
amis_mod <- inlaAMIS(data = df, init = init, prior.x.mis,
                               dq.x.mis, rq.x.mis, fit.inla,
                              N_t = seq(25,50,1)*10, N_0 = 250,ncores= 10)
save(amis_mod, file = "./sims/missing/amis_missing.Rdata")

### IS
set.seed(1)
is_mod <- inlaIS(data = df, init = init, prior.x.mis,
                           dq.x.mis, rq.x.mis,fit.inla, N_0 = 800, N = 10000,ncores = 10)
save(is_mod, file = "./sims/missing/is_missing.Rdata")

### MCMC
set.seed(1)
# proposal distribution
## evaluate
dq.x.mis <- function(y, eta, sigma=sqrt(5), log =TRUE) {
  sum(dnorm(y, mean = eta, sd = sigma, log = log))
}
## sample
rq.x.mis <- function(eta,sigma=sqrt(5)) {
  rnorm(length(eta), mean = eta, sd = sigma)
}
init = list(rep(mean(df$d.mis$bmi, na.rm = TRUE),n.mis), sqrt(5))
mcmc_mod <- inlaMH(data = df, init = init,
                               prior.x.mis, dq.x.mis, rq.x.mis, fit.inla,
                               n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_mod, file = "./sims/missing/mcmc_missing.Rdata")


