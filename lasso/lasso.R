# loading required packages
library(INLA)
library(ISLR)
library(glmnet)
library(INLA)
library(smoothmest)
library(mvtnorm)

# sourcing INLA-IS, INLA-AMIS and INLA-MH code
source("./inlaMC/inlaMC.R")

data(Hitters)

#Check NA's and fix

Hitters <- na.omit(Hitters)

#
# The Lasso
#

#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)
df <- list(y = y, x = x)
n.beta <- ncol(df$x)

# ml estimates
ml = summary(lm(y~-1 + x, data = df))$coefficients[,1:2]

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)

#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit lasso model for several values of lambda
lasso.mod <- glmnet(x[train, ] , y[train], alpha = 1, lambda = grid,intercept = F)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1,intercept=F)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)


#Fitted values
lasso.fitted <- predict(out, s = bestlam, newx = x)
# importing dataset

# finding inverse of the precision
stdev.samp <- .25 * solve(t(x)%*%x)


fit.inla <- function(data, eta) {
  data$oset = data$x %*% eta
  res = inla(y ~ -1 + offset(oset), data = data)
  res = inla.rerun(res)
  return(list(mlik = res$mlik[[1]],
              dists = list(tau = res$marginals.hyperpar[[1]]),
              stats = list(tau = as.numeric(res$summary.hyperpar[1]))))
}



prior.beta <- function(x, mu = 0, lambda = 0.073, log = TRUE) {
  res <- sum(log(ddoublex(x, mu = mu, lambda = lambda)))

  if(!log) { res <- exp(res) }

  return(res)
}

# initial parameters of the proposal distribution
init = list(mu = rep(0,n.beta), cov = 4*stdev.samp)
# proposal distribution
## evaluate
dq.beta <- function(y, theta = init, log =TRUE) {
  #dmvnorm(y,mean = x, sigma = sigma,log = log)
  dmvt(y,delta=theta[[1]],sigma=theta[[2]],df=3,log=log,type = "shifted")
}
## sample
rq.beta <- function(theta) {
  #rmvnorm(1,mean=x,sigma = sigma)
  as.vector(rmvt(1,sigma = theta[[2]], df=3, delta = theta[[1]], type = "shifted"))
}


### AMIS
set.seed(1)
amis_mod = inlaAMIS(data = df, init = init, prior.beta, dq.beta, 
                    rq.beta, fit.inla, N_t = seq(25,50,1)*10, N_0 = 25,ncores = 10)
save(amis_mod, file = "./sims/lasso/amis_lasso.Rdata")

### IS
set.seed(1)
is_mod = inlaIS(data = df, init = list(mu = rep(0,n.beta), cov = 4*stdev.samp),
                          prior.beta, dq.beta, rq.beta, fit.inla, N_0 = 800, N = 10000,ncores = 10)
save(is_mod, file = "./sims/lasso/is_lasso.Rdata")


### MCMC
dq.beta <- function(x, y, sigma = stdev.samp, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
  #rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted")
}

rq.beta <- function(x, sigma = stdev.samp) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
  #rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted")
}
set.seed(1)
mcmc_mod = inlaMH(data = df, init = list(mu = rep(x=0,ncol(df$x)),cov = stdev.samp),
                              prior.beta, dq.beta, rq.beta, fit.inla,
                              n.samples = 10500, n.burnin = 500, n.thin = 1)
save(mcmc_mod, file = "./sims/lasso/mcmc_lasso.Rdata")

