# loading required packages
library(INLA)
library(ISLR)
library(glmnet)
library(smoothmest)
library(mvtnorm)


# loading general functions
source("./genFuncs.R")
# loading amis with inla functions
source("./inlaAMIS.R")
# loading is with inla functions
source("./inlaIS.R")
# loading mcmc with inla functions
source("./inlaMH.R")

data(Hitters)
summary(Hitters)

#Check NA's and fix
sum(is.na(Hitters$Salary))
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
plot(lasso.mod)

#CV
set.seed(1)
cv.out <- cv.glmnet(x[train, ], y[train], alpha = 1,intercept=F)
plot(cv.out)

#Take best lambda for lasso model
bestlam <- cv.out$lambda.min

#Predcit with lasso on test data
lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test, ])
mean((lasso.pred - y[test])^2)

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = bestlam)

#Check estimated coefficients
lasso.coef
lasso.coef[lasso.coef != 0]

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

dq.beta <- function(y, x, sigma = stdev.samp, log =TRUE) {
  #dmvnorm(y,mean = x, sigma = sigma,log = log)
  dmvt(y,delta=x,sigma=sigma,df=3,log=log,type = "shifted")
}

rq.beta <- function(x, sigma = stdev.samp) {
  #rmvnorm(1,mean=x,sigma = sigma)
  as.vector(rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted"))
}


### AMIS
amis_mod = inlaAMIS(data = df, init = list(mu = rep(0,n.beta), cov = 4*stdev.samp), 
                    prior.beta, dq.beta, rq.beta, fit.inla, N_t = seq(25,50,1), N_0 = 25,ncores = 10)
save(amis_mod, file = "./lasso/lasso-amis-w-inla.Rdata")

### IS
is_mod = inlaIS(data = df, init = list(mu = rep(0,n.beta), cov = 4*stdev.samp), 
                          prior.beta, dq.beta, rq.beta, fit.inla, N_0 = 800, N = 10000,ncores = 10)
save(is_mod, file = "./lasso/lasso-is.Rdata")



### MCMC
source("./lasso/lasso_mcmc_w_inla.R")
dq.beta <- function(x, y, sigma = stdev.samp, log =TRUE) {
  dmvnorm(y, mean = x, sigma = sigma, log = log)
  #rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted")
}

rq.beta <- function(x, sigma = stdev.samp) {
  as.vector(rmvnorm(1, mean = x, sigma = sigma))
  #rmvt(1,sigma = sigma, df=3, delta = x, type = "shifted")
}
mcmc_mod = inlaMH(data = df, init = rep(0,ncol(df$x)), 
                              prior.beta, dq.beta, rq.beta, fit.inla, 
                              n.samples = 100500, n.burnin = 500, n.thin = 10)
save(mcmc_mod, file = "./lasso/lasso-mcmc-w-inla.Rdata")

