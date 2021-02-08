# calculating parameters of proposal distribution AMIS and IS
calc.theta <- function(theta,weight,eta,i_tot,i_cur){
  weight[1:i_tot] = exp(weight[1:i_tot] - max(weight[1:i_tot]))
  for (i in seq(ncol(eta))){
    theta$a.mu[i_cur,i] = sum(eta[1:i_tot,i]*weight[1:i_tot])/sum(weight[1:i_tot])
  }
  for (i in seq(ncol(eta))){
    for (j in seq(i,ncol(eta))){
      theta$a.cov[i,j,i_cur] = theta$a.cov[j,i,i_cur] = sum(weight[1:i_tot]*(eta[1:i_tot,i]-theta$a.mu[i_cur,i])*
                                                              (eta[1:i_tot,j]-theta$a.mu[i_cur,j]))/(sum(weight[1:i_tot]))
    }
  }
  return(theta)
}

# calculating mean and variance of posterior marginals AMIS and IS
calc.stats <- function(stats,weight){
  for (i in seq(length(stats))){
    new.stat = c(0,0)
    new.stat[1] =  sum(stats[[i]]*weight)/sum(weight)
    new.stat[2] = sum((stats[[i]] - new.stat[1])*(stats[[i]] - new.stat[1])*weight)/sum(weight)
    stats[[i]] = new.stat
  }
  return(stats)
}

# adding stat to a list
store.stats <- function(stat,stats,j,n.prop){
  if (anyNA(stats)){
    stats = stat
    for (i in seq(length(stat))){
      stats[[i]] = rep(NA, n.prop)
      stats[[i]][j] = stat[[i]]
    }
    return(stats)
  }else{
    for (i in seq(length(stat))){
      stats[[i]][j] = stat[[i]]
    }
    return(stats)
  }
}

# adding conditional posterior marginal densities to a list
store.post <- function(marg,margs,j,n.prop){
  if (anyNA(margs)){
    margs = marg
    for (i in seq(length(marg))){
      margs[[i]] = list(x = matrix(NA, nrow = n.prop, ncol = length(marg[[i]][,1])),
                        y = matrix(NA, nrow = n.prop, ncol = length(marg[[i]][,2])))
      margs[[i]]$x[j,] = marg[[i]][,1]
      margs[[i]]$y[j,] = marg[[i]][,2]
    }
    return(margs)
  }else{
    for (i in seq(length(marg))){
      margs[[i]]$x[j,] = marg[[i]][,1]
      margs[[i]]$y[j,] = marg[[i]][,2]
    }
    return(margs)
  }
}

kde_mc <- function(eta,weight){
  weight = weight/sum(weight)
  usable = weight>1e-10
  weight = weight[usable]
  eta = matrix(eta[usable,],ncol=ncol(eta))
  return(lapply(seq(ncol(eta)), function(x){
    as.data.frame(density(x = eta[,x],
                          weights = weight,
                          kernel = "gaussian")[c(1,2)])
  }))
}

# calculating the quantile curves using the true parameter values
pqr_truth_lines <- function(data,params,type){
  quants = c(0.1,0.25,0.5,0.75,0.9)
  mu = params[[1]] + params[[2]]*data$x
  tau = exp(params[[3]] + params[[4]]*data$x)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (i in seq(length(quants))){
    if (type == "gaussian"){
      tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
    }else if (type == "gamma"){
      tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau
    }
    res = rbind(res,data.frame(x = data$x, y = tmpy, quants = rep(toString(quants[i]),length(data$x))))
  }
  return(res[-1,])
}

# calculating the quantile curves using the approximated parameters
pqr_inla <- function(data,margs,a.mu,type,domain = NA){
  tmpx = data$x
  quants = c(0.1,0.25,0.5,0.75,0.9)
  params = c(INLA::inla.zmarginal(margs$a,silent=T)[[1]],
             INLA::inla.zmarginal(margs$b,silent=T)[[1]],
             a.mu[1],
             a.mu[2])
  if (!anyNA(domain)){
    tmpx = seq(domain[1],domain[2],length.out = 500)
  }
  mu = params[1] + params[2]*tmpx
  tau = exp(params[3] + params[4]*tmpx)
  res = data.frame(x = NA, y = NA, quants = NA)
  for (i in seq(length(quants))){
    if (type == "gaussian"){
      tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
    }else if (type == "gamma"){
      tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau
    }
    res = rbind(res,data.frame(x = tmpx, y = tmpy, quants = rep(toString(quants[i]),length(tmpx))))
  }
  return(res[-1,])
}

# function calculating poster mean
calc.post.mean <- function(eta,ws=NA){
  if (anyNA(ws)){
    colMeans(eta)
  }else{
    colSums(eta*ws)/sum(ws)
  }
}

# function calculating posterior standard deviation
calc.post.sd <- function(eta,ws=NA){
  res = numeric(ncol(eta))
  if (anyNA(ws)){
    for (i in seq(ncol(eta))){
      res[i] = sd(eta[,i])
    }
  }else{
    tmp = calc.post.mean(eta,ws)
    for (i in seq(ncol(eta))){
      res[i] = sqrt(sum(ws*(eta[,i] - tmp[i])^2/sum(ws)))
    }
  }
  return(res)
}

# calculating running effective sample size
running.ESS <- function(eta, times, ws = NA, norm = TRUE,step = 100){
  if (anyNA(ws)){
    require(coda)
    ess = sapply(lapply(seq(2,nrow(eta)),function(x){
      coda::effectiveSize(eta[1:x,])
    }),min)
    times = times[-1]
  }else{
    if (norm){
      ws = ws/sum(ws)
    }
    ess = sapply(seq(length(ws)),function(x){
      sum(ws[1:x])^2/(sum(ws[1:x]^2))
    })
    rm.ess = !is.na(ess)
    times = times[rm.ess]
    ess = ess[rm.ess]
  }
  ess.df = data.frame(time = c(times[1],times[rev(seq(length(times),100,-step))]),
                      ess = c(ess[1],ess[rev(seq(length(ess),100,-step))]))
  return(ess.df)
}

# Bayesian model averaging of conditional posterior marginals
fit.marginals <- function(ws,margs,len = 400){
  ws = ws/sum(ws)
  xmin <- quantile(apply(margs[[1]],1,function(X){min(X)}),0.25)
  xmax <- quantile(apply(margs[[1]],1,function(X){max(X)}),0.75)
  xx <- seq(xmin, xmax, len = len)
  marg = numeric(len)
  for (i in seq(nrow(margs[[1]]))){
    marg = marg + ws[i]*INLA::inla.dmarginal(xx, list(x = margs[[1]][i,], y = margs[[2]][i,]))
  }
  data.frame(x = xx, y = marg)
}

# multivariate weighted kernel density estimation
kde2d.weighted <- function (x, y, w, h, n = 25, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx)
    stop("data vectors must be the same length")
  if (length(w) != nx & length(w) != 1)
    stop("weight vectors must be 1 or length of data")
  gx <- seq(lims[1], lims[2], length = n)
  gy <- seq(lims[3], lims[4], length = n)
  if (missing(h))
    h <- c(MASS::bandwidth.nrd(x), MASS::bandwidth.nrd(y));
  if (missing(w))
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1]
  ay <- outer(gy, y, "-")/h[2]
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2])
  return(list(x = gx, y = gy, z = z))
}

