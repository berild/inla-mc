# main MCMC with INLA algorithm
inlaMH <- function(data, init, prior, d.prop, r.prop, fit.inla,
                    n.samples = 100, n.burnin = 5, n.thin = 1,kde = T){
  eta = matrix(data = NA,nrow = n.samples, ncol = length(init$mu))
  mlik = numeric(n.samples)
  acc.vec = numeric(n.samples)
  eta[1,] = init$mu
  mod.curr = fit.inla(data=data, eta = eta[1,])
  mlik[1] = mod.curr$mlik
  starttime = Sys.time()
  pb <- txtProgressBar(min = 0, max = n.samples, style = 3)
  i_marg = 0
  N_marg = floor((n.samples - n.burnin)/n.thin)
  times = numeric(N_marg)
  margs = NA
  for (i in seq(2,n.samples)){
    setTxtProgressBar(pb, i)
    INLA_crash = T
    while(INLA_crash){
      tryCatch({
        eta.new = r.prop(eta = eta[i-1,], sigma = init$cov)
        mod.new = fit.inla(data = data,eta = eta.new)
        INLA_crash = F
      },error=function(e){
      },finally={})
    }
    lacc1 = mod.new$mlik + prior(eta.new) + d.prop(eta.new, eta[i-1,],init$cov)
    lacc2 = mod.curr$mlik + prior(eta[i-1,]) + d.prop(eta[i-1,], eta.new,init$cov)
    acc = min(1,exp(lacc1 - lacc2))
    if (runif(1) < acc){
      eta[i,] = eta.new
      mod.curr = mod.new
      mlik[i] = mod.new$mlik
      acc.vec[i] = T
    }else{
      eta[i,] = eta[i-1,]
      mlik[i] = mlik[i-1]
      acc.vec[i] = F
    }
    if(i > n.burnin){
      if (((i-1) %% n.thin)==0){
        i_marg = i_marg + 1
        margs = store.post(mod.curr$dists,margs,i_marg,N_marg)
        times[i_marg] =  as.numeric(difftime(Sys.time(),starttime,units = "secs"))
      }
    }
  }
  eta = eta[-seq(n.burnin),]
  eta = eta[seq(from = 1, to = nrow(eta), by=n.thin),]
  res$eta = eta
  res$times = times
  res$mlik = mlik
  res$acc.vec = acc.vec
  res$margs = lapply(margs, function(x){fit.marginals(rep(1,nrow(eta)),x)})
  if (kde){
    kde_mc(res$eta,rep(1,nrow(eta)))
  }
  return(res)
}
