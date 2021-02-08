# function available for parallellization
par.is <- function(x, data, theta, t, prior, d.prop, r.prop, fit.inla){
  INLA_crash = T
  while(INLA_crash){
    tryCatch({
      eta = r.prop(theta$a.mu[t,], sigma = theta$a.cov[,,t])
      mod = fit.inla(data = data ,eta = eta)
      INLA_crash = F
    },error=function(e){
    },finally={})
  }
  weight = mod$mlik + prior(eta) - d.prop(y = eta, x = theta$a.mu[t,], sigma = theta$a.cov[,,t])
  return(list(mlik = mod$mlik, dists = mod$dists, eta = eta, weight = weight, times = Sys.time()))
}

# main IS with INLA algorithm
inlaIS <- function(data, init, prior, d.prop, r.prop, fit.inla, N_0 = NA, N = 4000, kde = TRUE, ncores=1){
  if (ncores>parallel::detectCores()){
    ncores = parallel::detectCores()
  }
  theta = list(a.mu = matrix(NA, ncol = length(init$mu), nrow = 2),
               a.cov = array(NA, dim = c(length(init$mu), length(init$mu), 2)))
  theta$a.mu[1,] = init$mu
  theta$a.cov[,,1] = init$cov
  starttime = Sys.time()
  times = numeric(N)
  res = list()
  if (anyNA(N_0)){
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    theta$a.mu[2,] = init$mu
    theta$a.cov[,,2] = init$cov
    N_0 = 0
  }else{
    pb <- txtProgressBar(min = 0, max = N+N_0, style = 3)
    weight = numeric(N_0)
    eta = matrix(NA, ncol = length(init$mu), nrow = N_0)
    is.list = parallel::mclapply(seq(N_0), function(x){
      par.is(x, data, theta, 1, prior, d.prop,r.prop, fit.inla)
    }, mc.set.seed = TRUE, mc.cores = ncores)
    for (i in seq(length(is.list))){
      setTxtProgressBar(pb, i)
      eta[i,]= is.list[[i]]$eta
      weight[i] = is.list[[i]]$weight
    }
    theta = calc.theta(theta,weight,eta,N_0,2)
  }
  margs = NA
  eta = matrix(NA, ncol = length(init$mu), nrow = N)
  weight = numeric(N)
  mlik = numeric(N)
  is.list = mclapply(seq(N), function(x){
    par.is(x, data, theta, 2, prior, d.prop,r.prop, fit.inla)
  }, mc.set.seed = TRUE, mc.cores = ncores)
  for (i in seq(length(is.list))){
    setTxtProgressBar(pb, i+N_0)
    eta[i,] = is.list[[i]]$eta
    margs = store.post(is.list[[i]]$dists,margs,i,N)
    weight[i] = is.list[[i]]$weight
    mlik[i] = is.list[[i]]$mlik
    times[i] = as.numeric(difftime(is.list[[i]]$times,starttime,units = "secs"))
  }
  res$eta = eta
  res$times = times
  res$theta = theta
  res$weight = exp(weight - max(weight))
  res$margs = lapply(margs, function(x){fit.marginals(res$weight,x)})
  if (kde){
    res$eta_kern = kde_mc(res$eta,res$weight)
  }
  return(res)
}
