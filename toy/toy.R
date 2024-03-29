library(mvtnorm)
library(MASS)
library(INLA)

# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("./inlaMC/inlaMC.R")

# function to sample data
sample.linreg <- function(){
  n = 100
  x1 = runif(n)
  x2 = runif(n)
  err = rnorm(n)
  y = 1 + 1*x1 -1*x2 + err
  return(list(y = y,x = matrix(c(x1,x2),ncol = 2)))
}

# sampling dataset
set.seed(1)
df = sample.linreg()

# fit the exact inla model
inla_mod = inla(y~x, data=df) 
save(inla_mod, file = "./sims/toy/inla-toy.Rdata")

# inla function fitting conditional LGMs in the combined methods
fit.inla <- function(data, eta){
  data$oset = data$x%*%eta
  res = inla(y~1+offset(oset), data = data)
  return(list(mlik = res$mlik[1],
              dists = list(intercept = res$marginals.fixed[[1]], 
                           tau = res$marginals.hyperpar[[1]])))
}

# the prior for the beta parameters
prior.beta <- function(x, sigma = sqrt(1/.001), log = TRUE) {
  sum(dnorm(x, mean = 0, sd= sigma, log = log))
}

# initial parameters of the proposal distribution in AMIS and IS
init = list(mu = c(0,0),cov = diag(5,2,2))

# proposal distribution in AMIS and IS
rq.beta <- function(theta = init) {
  as.vector(rmvnorm(1,mean=theta[[1]],sigma = theta[[2]]))
}

dq.beta <- function(y, theta = init, t = 1, log =TRUE) {
  dmvnorm(y,mean = theta[[1]], sigma = theta[[2]],log = log)
}

### AMIS
set.seed(1)
amis_mod <- inlaAMIS(data = df, init = init, prior.beta,
                               dq.beta, rq.beta, fit.inla,
                               N_t = seq(25,50)*1, N_0 = 2,ncores = 10,kde = T)
save(amis_mod, file = "./sims/toy/amis_toy.Rdata")

### IS
set.seed(1)
is_mod <- inlaIS(data = df, init = init, prior.beta,
                           dq.beta, rq.beta, fit.inla, N_0 = 800, N = 10000,ncores = 10,kde = T)
save(is_mod, file = "./sims/toy/is_toy.Rdata")

### MCMC

# proposal distribution MCMC
rq.beta <- function(x, sigma = .75) {
  rnorm(length(x), mean = x, sd = sigma)
}
dq.beta <- function(y, x, sigma = .75, log =TRUE) {
  sum(dnorm(x, mean = y, sd = sigma, log = log))
}

# initial state MCMC
init = list(c(0,0), .75)
set.seed(1)
mcmc_mod <- inlaMH(data = df, init = init, prior.beta,
                               dq.beta, rq.beta, fit.inla, n.samples = 10500, n.burnin = 500, n.thin = 1, kde=T)
save(mcmc_mod, file = "./sims/toy/mcmc_toy.Rdata")
