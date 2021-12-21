# function available for parallellization
par.is <- function(x, data, theta, t, prior, d.prop, r.prop, fit.inla){
  INLA_crash = T
  while(INLA_crash){
    tryCatch({
      if (t==1){
        eta = r.prop(theta)
      }else{
        eta = r.prop(theta[[2]]) 
      }
      mod = fit.inla(data ,eta)
      INLA_crash = F
    },error=function(e){
    },finally={})
  }
  if (t==1){
    weight = mod$mlik + prior(eta) - d.prop(eta, theta)
  }else{
    weight = mod$mlik + prior(eta) - d.prop(eta, theta[[2]])
  }
  return(list(mlik = mod$mlik, dists = mod$dists, eta = eta, weight = weight, times = Sys.time()))
}

# main IS with INLA algorithm
inlaIS <- function(data, init, prior, d.prop, r.prop, fit.inla, N_0 = NA, N = 4000, kde = FALSE, ncores=1){
  if (ncores>parallel::detectCores()){
    ncores = parallel::detectCores()
  }
  theta = init
  starttime = Sys.time()
  times = numeric(N)
  res = list()
  margs = NA
  n_eta = length(r.prop(theta))
  if (anyNA(N_0)){
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    theta = list(theta, init)
    N_0 = 0
  }else{
    pb <- txtProgressBar(min = 0, max = N+N_0, style = 3)
    weight = numeric(N_0)
    eta = matrix(NA, ncol = n_eta, nrow = N_0)
    is.list = parallel::mclapply(seq(N_0), function(x){
      par.is(x, data, theta, 1, prior, d.prop,r.prop, fit.inla)
    }, mc.set.seed = TRUE, mc.cores = ncores)
    for (i in seq(length(is.list))){
      setTxtProgressBar(pb, i)
      eta[i,]= is.list[[i]]$eta
      weight[i] = is.list[[i]]$weight
      margs = store.post(is.list[[i]]$dists,margs,i,N)
    }
    theta = list(theta,calc.theta(list(weight,eta,margs),N_0))
  }
  margs = NA
  
  eta = matrix(NA, ncol = n_eta, nrow = N)
  weight = numeric(N)
  mlik = numeric(N)
  is.list = parallel::mclapply(seq(N), function(x){
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
  res$dists = res$margs
  res$margs = lapply(margs, function(x){fit.marginals(res$weight,x)})
  if (kde){
    res$eta_kern = kde_mc(res$eta,res$weight)
  }
  return(res)
}
