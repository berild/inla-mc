library(parallel)
library(Brq)
library(INLA)
library(SemiPar)

# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("./inlaMC/inlaMC.R")

# function fitting conditional models with INLA with gamma likelihood
fit.inla.ggamma <- function(data,eta){
  res = inla(y ~ 1 + x,
             data = data,
             scale = exp(eta[1] + eta[2]*data$x),
             family = "gamma",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(a = res$marginals.fixed[[1]],
                           b = res$marginals.fixed[[2]])))
}

# function fitting conditional models with INLA with gaussian likelihood
fit.inla.gaussian <- function(data,eta){
  res = inla(y ~ 1 + x,
             data = data,
             family = "gaussian",
             scale = exp(eta[1] + eta[2]*data$x),
             control.family = list(hyper = list(prec = list(initial = log(1), fixed = TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(a = res$marginals.fixed[[1]],
                           b = res$marginals.fixed[[2]])))
}

# function fitting conditional RW2 models with INLA
fit.inla.rw2 <- function(data,eta){
  formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
  res = inla(formula,
             data = data,
             control.predictor = list(compute = T),
             scale = exp(eta[1] + eta[2]*data$range),
             family = "gaussian",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           prec.range = res$marginals.hyperpar[[1]])))
}

# function fitting conditional RW2 models with INLA
fit.inla.rw2.gamma <- function(data,eta){
  formula = IgG ~ f(Age, model = "rw2", constr = T, scale.model =  F)
  res = inla(formula,
             data = data,
             control.predictor = list(compute = T),
             scale = exp(eta[1] + eta[2]*data$Age),
             family = "gamma",
             control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
             verbose = FALSE)
  return(list(mlik = res$mlik[[1]],
              dists = list(intercept = res$marginals.fixed[[1]],
                           prec.range = res$marginals.hyperpar[[1]])))
}

# vague prior for the parameters
prior.param <- function(x, log = TRUE) {
  sum(dnorm(x, 0, 100, log = log))
}

# changing intital proposal distribution amis and is
init = list(mu = c(10,-10),cov = 3*diag(2))

# proposal distribution amis and is
dq.param <- function(y, theta = init, log =TRUE) {
  mvtnorm::dmvt(y,sigma = theta[[2]], df=3, delta = theta[[1]], type = "shifted",log=log)
  #dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.param <- function(theta=init) {
  as.vector(mvtnorm::rmvt(1,sigma = theta[[2]], df=3, delta = theta[[1]], type = "shifted"))
  #as.vector(rmvnorm(1, mean = x, sigma = sigma))
}

# Lidar data
data(lidar)
#scaling data for easier convergence with current proposal
r_max = max(lidar$range)
r_min = min(lidar$range)
lidar$range = lidar$range/r_max


# running amis with inla on lidar data
set.seed(1)
amis_mod = inlaAMIS(data = lidar, init = init, prior.param,
                              dq.param, rq.param, fit.inla.rw2,
                              N_t = seq(25,50,1)*10, N_0 = 250,ncores = 10)
amis_w_inla_mod$mod = lidar
amis_w_inla_mod$scale = c(r_max,r_min)
save(amis_w_inla_mod, file = "./sims/pqr/amis_pqr_rw2.Rdata")

# running is with inla on lidar data
set.seed(1)
is_mod = inlaIS(data = lidar, init = init, prior.param,
                          dq.param, rq.param, fit.inla.rw2,
                          N_0 = 800, N = 10000,ncores = 10)
is_w_inla_mod$mod = lidar
is_w_inla_mod$scale = c(r_max,r_min)
save(is_w_inla_mod, file = "./sims/pqr/is_pqr_rw2.Rdata")

# initalizing proposal distribution for mcmc with inla
init_mcmc = list(mu = c(0,0),cov = diag(2))
dq.param.mcmc <- function(y, x, sigma = init$cov, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
}

rq.param.mcmc  <- function(x, sigma = init$cov) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
}
# running mcmc with inla on lidar data
set.seed(1)
mcmc_mod <- inlaMH(data = lidar, init = init_mcmc,
                               prior.param, dq.param.mcmc, rq.param.mcmc, fit.inla.rw2,
                               n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_w_inla_mod, file = "./sims/pqr/mcmc_pqr_rw2.Rdata")
eta_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1,nrow(mcmc_w_inla_mod$eta))/nrow(mcmc_w_inla_mod$eta), n = 100, lims = c(-1,1,-1,1))
mcmc_w_inla_mod$eta_kern_joint = data.frame(expand.grid(x=eta_kern_mcmc$x, y=eta_kern_mcmc$y), z=as.vector(eta_kern_mcmc$z))
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                        kernel = "gaussian")[c(1,2)])
})
mcmc_w_inla_mod$mod = lidar
mcmc_w_inla_mod$scale = c(r_max,r_min)
save(mcmc_w_inla_mod, file = "./sims/pqr-gaussian-lidar-mcmc-w-inla.Rdata")

# ImmunogG data
data("ImmunogG")
# changing intial proposal distribution amis and is

init = list(mu = c(2,0),cov = diag(c(1,0.5)))

formula = IgG ~ f(Age, model = "rw2", constr = T, scale.model =  F)
res = inla(formula,
           data = ImmunogG,
           control.predictor = list(compute = T),
           family = "gamma",
           verbose = FALSE)

#running amis with inla on ImmunogG data
amis_mod = inlaAMIS(data = ImmunogG,
                    init = init, prior.param, dq.param, rq.param, fit.inla.rw2.gamma,
                    N_t = seq(25,50,1), N_0 = 25,kde = T,ncores=10)
save(amis_mod, file = "./pqr/pqr-gammarw2-ImmunogG-amis.Rdata")


# running is with inla on ImmunogG data
is_w_inla_mod = is.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG),
                         init = init, prior.param, dq.param, rq.param, fit.inla.ggamma, N_0 = 800,N = 10000)
is_w_inla_mod$mod = list(x=ImmunogG$Age,y = ImmunogG$IgG)
save(is_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-is-w-inla.Rdata")

# running mcmc with inla on ImmunogG data
init_mcmc = list(mu = c(0,0),cov = diag(c(0.3,0.15)))
mcmc_w_inla_mod <- mcmc.w.inla(data = list(x=ImmunogG$Age,y = ImmunogG$IgG), init = init_mcmc,
                              prior.param, dq.param.mcmc, rq.param.mcmc, fit.inla.ggamma,
                              n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-mcmc-w-inla.Rdata")
eta_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1,nrow(mcmc_w_inla_mod$eta))/nrow(mcmc_w_inla_mod$eta), n = 100, lims = c(-1,1,-1,1))
mcmc_w_inla_mod$eta_kern_joint = data.frame(expand.grid(x=eta_kern_mcmc$x, y=eta_kern_mcmc$y), z=as.vector(eta_kern_mcmc$z))
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
 as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                       weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                       kernel = "gaussian")[c(1,2)])
})
mcmc_w_inla_mod$mod = list(x=ImmunogG$Age,y = ImmunogG$IgG)
save(mcmc_w_inla_mod, file = "./sims/pqr-gamma-ImmunogG-mcmc-w-inla.Rdata")
