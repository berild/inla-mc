# calc help parameter in weight calculation
calc.delta <- function(N_t,eta,theta,t,d.prop){
  tmp = 0
  for (l in seq(t)){
    #tmp = tmp + N_t[l]*d.prop(y = eta, x = theta$a.mu[l+1,], sigma = theta$a.cov[,,l+1], log = FALSE)
    tmp = tmp + N_t[l]*d.prop(eta, theta[[l+1]], log = FALSE)
  }
  return(tmp)
}

# function updating help parameter and weights
update.delta.weight <- function(delta,weight,N_t,eta,theta,t,mlik,prior,d.prop){
  i_tmp = 0
  N_tmp = sum(N_t[1:(t+1)])
  for (l in seq(t)){
    for (i in seq(N_t[l])){
      i_tmp = i_tmp + 1
      #delta[i_tmp] = delta[i_tmp] + N_t[l]*d.prop(y = eta[i_tmp,], x = theta$a.mu[t+1,], sigma = theta$a.cov[,,t+1], log = FALSE)
      delta[i_tmp] = delta[i_tmp] + N_t[l]*d.prop(eta[i_tmp,], theta[[t+1]], log = FALSE)
      weight[i_tmp] = mlik[i_tmp] + prior(eta[i_tmp,]) - log(delta[i_tmp]/N_tmp)
    }
  }
  return(list(
    delta = delta,
    weight = weight
  ))
}

# function being parallelized
par.amis <- function(x,data, theta, t, N_0, N_t, N_tmp,
                     prior, d.prop, r.prop, fit.inla){
  INLA_crash = TRUE
  while(INLA_crash){
    tryCatch({
      eta = r.prop(theta[[t+1]]) 
      mod = fit.inla(data ,eta)
      INLA_crash = F
    },error=function(e){
      system("echo 'Encoutered error in r.prop or fit.inla'")
    },finally={})
  }
  if (t==0){
    delta = N_0*d.prop(eta, theta[[1]], log = FALSE)
    weight = mod$mlik + prior(eta) - d.prop(eta, theta[[1]])
  }else{
    #delta = N_0*d.prop(y = eta, x = theta$a.mu[1,], sigma = theta$a.cov[,,1],log = FALSE) + calc.delta(N_t,eta,theta, t, d.prop)
    delta = N_0*d.prop(eta, theta[[1]], log = FALSE) + calc.delta(N_t,eta,theta, t, d.prop)
    weight = mod$mlik + prior(eta)- log(delta/N_tmp)
  }
  return(list(mlik = mod$mlik, dists = mod$dists, eta = eta, delta = delta, weight = weight, times = Sys.time()))
}



# main AMIS w/ INLA algorithm
inlaAMIS <- function(data, init, prior, d.prop, r.prop, fit.inla, N_t = seq(25,50)*10, N_0 = 250, kde = FALSE, ncores = 1){
  if (anyNA(N_0)){
    N_0 = round(sum(N_t)/2)
  }
  if (ncores>parallel::detectCores()){
    ncores = parallel::detectCores()
  }
  N_tot = N_0 + sum(N_t)
  mlik = numeric(N_tot)
  delta = numeric(N_tot)
  weight = numeric(N_tot)
  times = numeric(N_tot)
  theta  = list()
  theta[[1]] = init
  n_eta = length(r.prop(theta[[1]]))
  eta = matrix(NA, ncol = n_eta, nrow = N_tot)
  i_tot = 0
  pb <- txtProgressBar(min = 0, max = N_tot, style = 3)
  margs = NA
  starttime = Sys.time()
  N_tmp = N_0
  t = 0
  res = list()
  amis.list = parallel::mclapply(seq(N_0),function(x){
    par.amis(x, data, theta, t, N_0,
             N_t, N_tmp, prior, d.prop,
             r.prop, fit.inla)
  }, mc.cores = ncores)
  for (ele in amis.list){
    setTxtProgressBar(pb, i_tot)
    i_tot = i_tot + 1
    margs = store.post(ele$dists,margs,i_tot,N_tot)
    eta[i_tot,] = ele$eta
    mlik[i_tot] = ele$mlik
    delta[i_tot] = ele$delta
    weight[i_tot] = ele$weight
    times[i_tot] = as.numeric(difftime(ele$times,starttime,units = "secs"))

  }
  theta[[2]] = calc.theta(list(weight,eta,margs),i_tot)
  # adaptive importance sampling
  for (t in seq(length(N_t))){
    N_tmp = N_tmp + N_t[t]
    amis.list = parallel::mclapply(seq(N_t[t]),function(x){
      par.amis(x, data, theta, t, N_0,
               N_t, N_tmp, prior, d.prop,
               r.prop, fit.inla)
    }, mc.cores = ncores)
    for (ele in amis.list){
      setTxtProgressBar(pb, i_tot)
      i_tot = i_tot + 1
      margs = store.post(ele$dists,margs,i_tot,N_tot)
      eta[i_tot,] = ele$eta
      mlik[i_tot] = ele$mlik
      delta[i_tot] = ele$delta
      weight[i_tot] = ele$weight
      times[i_tot] = as.numeric(difftime(ele$times,starttime,units = "secs"))
    }
    delta.weight = update.delta.weight(delta[1:(N_tmp - N_t[t])],
                                       weight[1:(N_tmp - N_t[t])],
                                       N_t = c(N_0,N_t),
                                       matrix(eta[1:(N_tmp - N_t[t]),],ncol = n_eta),
                                       theta,
                                       t,
                                       mlik[1:(N_tmp - N_t[t])],
                                       prior,
                                       d.prop)
    delta[1:(N_tmp - N_t[t])] = delta.weight$delta
    weight[1:(N_tmp - N_t[t])] = delta.weight$weight
    theta[[t+2]] = calc.theta(list(weight,eta,margs),i_tot)
  }
  res$eta = eta
  res$times = times
  res$theta = theta
  res$weight = exp(weight - max(weight))
  res$mod=data
  res$dists = margs
  res$margs = lapply(margs, function(x){fit.marginals(res$weight,x)})
  if (kde){
    res$eta_kern = kde_mc(res$eta,res$weight)
  }
  return(res)
}

