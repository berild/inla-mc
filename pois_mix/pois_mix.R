library("MixtureInf")
library("INLA")
library("parallel")
library("ggplot2")
library("rjags")
# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("./inlaMC/inlaMC.R")

data(earthquake)

#Add year
earthquake$year <- 1900:2006

#Summary of the data
summary(earthquake)


# Values of indicator variable z: 1, 2

d <- data.frame(y = earthquake$number)
n <- nrow(d)

# Set seed
inits.jags <- list(
  z = sample(1:2, nrow(earthquake), rep = TRUE), # Initial groups
  ".RNG.name"="base::Mersenne-Twister", ".RNG.seed" = 1) # Set seed
  


# Fit model using JAGS (for checking)
jm <- jags.model("cluster.bug",
  data = list(y = earthquake$number, n = nrow(earthquake)),
  inits = inits.jags
)

# Burn-in
update(jm, n.iter = 10000)

# Sample from model
js <- jags.samples(jm,
  variable.names = c("z", "mu"),
  n.iter = 10000, n.thin = 10)

# Summary stats
apply(js$mu[ , ,1], 1, mean)
apply(js$mu[ , ,1], 1, sd)


# save results
save(file = "sims/pois-mix/mcmc_pois_mix.RData", list = c("jm", "js"))

#Summary 
apply(js$mu, 1, mean)
apply(js$mu, 1, sd)

# Proportion in group 1
summary(apply(js$z[,,1] == 1, 2, mean))
# Proportion in group 2
summary(apply(js$z[,,1] == 2, 2, mean))

# FIt model using IS/AMIS with INLA


fit.inla <- function(data, eta) {
  Y <- matrix(NA, ncol = 2, nrow = n)
  
  for(i in 1:2) {
    Y[which(eta == i), i] <- data$y[which(eta == i)]
  }
  
  # Create 2-column intercept
  Intercept <- Y
  Intercept[!is.na(Y)] <- 1
  
  res <- inla(Y ~ -1 + Intercept, data = list(Y = Y, Intercept = Intercept),
             family = rep("poisson", 2),
             control.fixed = list(mean = list(Intercept1 = log(10),
                                              Intercept2 = log(30)), prec = 0.01),
             num.threads = "1:1")
  return(list(mlik = res$mlik[[1]],
              dists = list(Intercept1 = res$marginals.fixed$Intercept1,
                           Intercept2 = res$marginals.fixed$Intercept2)))
}


ggplot(data.frame(x= d$y)) + 
  geom_histogram(aes(x=x),bins = 9,color="black",fill = "gray69")+ 
  labs(y="") + 
  theme_classic()
# Check

tmp <- fit.inla(d, sample(1:2, n, rep = TRUE))



# Get initial values of w and means using k-means
k_inits <- kmeans(d$y, 2)


# We assume that mean group 1 < mean group 2 (identifiability)
idx <- order(k_inits$centers)

# Initial proportions
w_init <- (k_inits$size / sum(k_inits$size)) [idx]
# We assume that mean group 1 < mean group 2 (identifiability)
means_init <- k_inits$centers[idx]

# initialize the poisson mean and mixture weight of proposal
init = list(w = w_init, mu = means_init)
# Sample values of z
# Assume two groups
# d: data; response in 'y'
# w: proportions of observations in each group
# x: Means of each group
rq.prop <- function(theta = init, y = d$y){
  w = theta$w
  x = theta$mu
  liks <- lapply(1:2, function(X) {
    w[X] * dpois(y, x[X])
  })
  
  # Probabilty of being in group 1
  probs <- liks[[1]] / (liks[[1]] + liks[[2]]) 
  
  z <- (runif(n) > probs) + 1
  
  return(z)
}

# Compute density using sampling distribution
# z: Vector of (sampled) indices
# y: data; response in 'y'
# w: proportions of observatiosn in each group
# x: Means of each group
dq.prop <- function(z, theta, log = TRUE, y = d$y) {
  # Compute probabilities
  x = theta$mu
  w = theta$w
  aux <- cbind(w[1] * dpois(y, x[1]), w[2] * dpois(y, x[2]))
  probs <- sapply(1:length(y), function(X) {
    (aux[X,] / sum(aux[X, ]))[z[X]]
  })
  
  if(log) {
    return(sum(log(probs)))
  }
  else {
    return(prod(probs))
  }
}

# Number of samples
n_samples <- 10000

prior <- function(eta){
  n*log(0.5)
}


# Adapting the mixture weight and poisson mean
calc.theta <- function(curr,i_tot){
  weight = curr[[1]]
  eta = curr[[2]]
  margs = curr[[3]]
  theta = list(mu = numeric(2),w = numeric(2)) 
  weight = weight[1:i_tot]
  weight = exp(weight - max(weight))
  weight = weight/sum(weight)
  marg1 = fit.marginals(weight,list(x = margs[[1]]$x[1:i_tot,],y = margs[[1]]$y[1:i_tot,]))
  marg2 = fit.marginals(weight,list(x = margs[[2]]$x[1:i_tot,],y = margs[[2]]$y[1:i_tot,]))
  theta$mu[1] = inla.emarginal(exp,marginal = marg1)
  theta$mu[2] = inla.emarginal(exp,marginal = marg2)
  theta$w[2] = sum(weight * (apply(eta[1:i_tot,], 1, mean) -1))
  theta$w[1] = 1 - theta$w[2]
  return(theta)
}


# IS-INLA call
set.seed(1)
is_mod <- inlaIS(d, init, prior = prior,d.prop = dq.prop, r.prop = rq.prop, 
                 fit.inla=fit.inla, N = 10, ncores = 60)
save(is_mod, file = "./sims/pois-mix/is_pois_mix.Rdata")

# AMIS-INLA call
set.seed(1)
amis_mod <- inlaAMIS(d, init=init, prior = prior, d.prop = dq.prop, r.prop = rq.prop,
                     fit.inla=fit.inla, N_t = rep(20, 4), N_0 = 10, ncores = 60)
save(amis_mod, file = "./sims/pois-mix/amis_pois_mix.Rdata")

sessionInfo()
