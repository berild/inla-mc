library(ggplot2)
library(INLA)
library(ggpubr)
source("./inlaMC/inlaMC.R")

# Loading results of simulations
load("./sims/pqr/mcmc_pqr_rw2.Rdata")
load("./sims/pqr/is_pqr_rw2.Rdata")
load("./sims/pqr/amis_pqr_rw2.Rdata")

# Rescaling the logratio fixed effect
# AMIS-INLA
domain = matrix(c(11,17,-0.022,-0.008),nrow=2)
amis_mod$theta[[1]][,2]=amis_mod$theta[[1]][,2]/amis_mod$scale[1]
amis_mod$params = amis_mod$theta[[1]][nrow(amis_mod$theta[[1]]),]
amis_mod$eta[,2] = amis_mod$eta[,2]/amis_mod$scale[1]
amis_mod$eta_kern = kde_mc(amis_mod$eta, amis_mod$weight)
amis_mod$mod$range = amis_mod$mod$range*amis_mod$scale[1]
# IS-INLA
is_mod$eta[,2] = is_mod$eta[,2]/is_mod$scale[1]
is_mod$params = calc.post.mean(is_mod$eta,is_mod$weight)
is_mod$eta_kern = kde_mc(is_mod$eta,is_mod$weight)
is_mod$mod$range = is_mod$mod$range*is_mod$scale[1]
# MCMC-INLA
mcmc_mod$eta[,2] = mcmc_mod$eta[,2]/mcmc_mod$scale[1]
mcmc_mod$params = calc.post.mean(mcmc_mod$eta)
mcmc_mod$eta_kern = kde_mc(mcmc_mod$eta, rep(1,nrow(mcmc_mod$eta)))
mcmc_mod$mod$range = amis_mod$mod$range*amis_mod$scale[1]

amis_mod$ess = running.ESS(amis_mod$eta, amis_mod$times,ws =  amis_mod$weight, norm = T)
is_mod$ess = running.ESS(is_mod$eta, is_mod$times,ws =  is_mod$weight, norm = T)
mcmc_mod$ess = running.ESS(mcmc_mod$eta, mcmc_mod$times)

# Function constructing Figure 3 in Supplementary material
figure3sup <- function(){
  height = 3
  width = 5
  p1 <- ggplot() +
    geom_line(data = as.data.frame(amis_mod$margs$intercept), aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data = as.data.frame(is_mod$margs$intercept),aes(x=x,y=y,linetype="IS with INLA")) +
    geom_line(data = as.data.frame(mcmc_mod$margs$intercept), aes(x=x,y=y,linetype="MCMC with INLA")) +
    labs(y = "",x=expression(mu[0]),linetype = "",title = "a") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.32,-0.26)) +
    theme_bw() +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))

  p3 <- ggplot() +
    geom_line(data = amis_mod$eta_kern[[1]],aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data = mcmc_mod$eta_kern[[1]],aes(x=x,y=y,linetype="MCMC with INLA")) +
    geom_line(data = is_mod$eta_kern[[1]],aes(x=x,y=y,linetype="IS with INLA")) +
    labs(y="",x = expression(alpha),linetype = "",title="b") +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    coord_cartesian(xlim=c(11.5,16)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
  p4 <- ggplot() +
    geom_line(data = amis_mod$eta_kern[[2]],aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data = mcmc_mod$eta_kern[[2]],aes(x=x,y=y,linetype="MCMC with INLA")) +
    geom_line(data = is_mod$eta_kern[[2]],aes(x=x,y=y,linetype="IS with INLA")) +
    labs(y= "",x = expression(beta),linetype = "",title = "c") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim=c(-0.018,-0.010)) +
    theme_bw() +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
  essp <- ggplot() +
    geom_line(data = amis_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA")) +
    geom_line(data = is_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA")) +
    geom_line(data = mcmc_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA")) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h","2h"),trans="log",breaks=c(0,60,60*5,60*20,60*60,2*60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title="d") +
    theme_bw() +
    coord_cartesian(xlim = c(10,2*60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12), axis.title = element_text(size=14))
  ggsave(filename = "pqr_uni_intercept.pdf", plot = p1, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "pqr_uni_alpha.pdf", plot = p3, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "pqr_uni_beta.pdf", plot = p4, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "pqr_uni_ess.pdf", plot = essp, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(list(p1,p3,p4,essp))
}
fig3sup <- figure3sup()
ggarrange(fig3sup[[1]],fig3sup[[2]],fig3sup[[3]],fig3sup[[4]],ncol=2, nrow=2, common.legend = T,legend="bottom")

# Calculating joint density of alpha and beta using kernel density estimation
# AMIS-INLA
eta_joint_kern_amis = kde2d.weighted(x = amis_mod$eta[,1], y = amis_mod$eta[,2], w = amis_mod$weight/(sum(amis_mod$weight)), n = 200, lims = c(11,16,-0.018,-0.009))
amis_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))
# IS-INLA
eta_joint_kern_is = kde2d.weighted(x = is_mod$eta[,1], y = is_mod$eta[,2], w = is_mod$weight/(sum(is_mod$weight)), n = 200, lims = c(11,16,-0.018,-0.009))
is_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_is$x, y=eta_joint_kern_is$y), z=as.vector(eta_joint_kern_is$z))
# MCMC-INLA
eta_joint_kern_mcmc = kde2d.weighted(x = mcmc_mod$eta[,1], y = mcmc_mod$eta[,2], w = rep(1/length(mcmc_mod$eta[,1]),length(mcmc_mod$eta[,1])), n = 200, lims = c(11,16,-0.018,-0.009))
mcmc_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_mcmc$x, y=eta_joint_kern_mcmc$y), z=as.vector(eta_joint_kern_mcmc$z))

# Function constructing Figure 4 in Supplementary material
figure4sup <- function(){
  nbins = 12
  height = 3
  width = 5
  p5 <- ggplot() +
    geom_contour(data = amis_mod$eta_joint_kern, aes(x = x, y = y, z = z,linetype = "AMIS with INLA"),bins = nbins, color = "black") +
    geom_contour(data = is_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,linetype = "IS with INLA"),bins = nbins,color = "black") +
    geom_contour(data = mcmc_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,linetype = "MCMC with INLA"),bins = nbins,color = "black") +
    theme_bw() +
    coord_cartesian(xlim = c(12.3,15.2), y = c(-0.0165,-0.0118)) +
    scale_linetype_manual(values = c(1,2,3))+
    labs(x=expression(alpha),y=expression(beta),linetype ="")+
    theme(legend.position="none",axis.title = element_text(size=12))
  p6 <- ggplot() +
    geom_contour(data = amis_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,linetype = "AMIS with INLA"),bins = nbins, color = "black") +
    geom_contour(data = is_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,linetype = "IS with INLA"),bins = nbins, color = "black") +
    geom_contour(data = mcmc_mod$eta_joint_kern, aes(x = x, y = y, z = z,linetype = "MCMC with INLA"),bins = nbins, color = "black") +
    theme_bw() +
    coord_cartesian(xlim = c(12.3,15.2), y = c(-0.0165,-0.0118)) +
    scale_linetype_manual(values = c(1,2,3))+
    labs(x=expression(alpha),y=expression(beta),linetype ="")+
    theme(legend.position="none",axis.title = element_text(size=12))
  p7 <- ggplot() +
    geom_contour(data = amis_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,linetype = "AMIS with INLA"),bins = nbins, color = "black") +
    geom_contour(data = mcmc_mod$eta_joint_kern, aes(x = x+60, y = y+60, z = z+60,linetype = "MCMC with INLA"),bins = nbins, color = "black") +
    geom_contour(data = is_mod$eta_joint_kern, aes(x = x, y = y, z = z,linetype = "IS with INLA"),bins = nbins, color = "black") +
    theme_bw() +
    coord_cartesian(xlim = c(12.3,15.2), y = c(-0.0165,-0.0118)) +
    scale_linetype_manual(values = c(1,2,3))+
    labs(x=expression(alpha),y=expression(beta),linetype ="")+
    theme(legend.position="none",axis.title = element_text(size=12))
  ggsave(filename = "pqr_joint_amis.pdf", plot = p5, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "pqr_joint_mcmc.pdf", plot = p6, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "pqr_joint_is.pdf", plot = p7, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(list(p5,p6,p7))
}
fig4sup <- figure4sup()
ggarrange(fig4sup[[1]],fig4sup[[2]],fig4sup[[3]],ncol=2, nrow=2, common.legend = T,legend="bottom")

# Getting Quantile Curves using INLA for AMIS-INLA fit
formula = logratio ~ f(range, model = "rw2", constr = T, scale.model =  F)
res_amis = inla(formula,
           data = amis_mod$mod,
           control.predictor = list(compute = T),
           scale = exp(amis_mod$params[1]+ amis_mod$params[2]*amis_mod$mod$range),
           control.family=list(hyper=list(theta=list(initial=log(1),fixed=TRUE))),
           verbose = FALSE)

figure5sup <- function(){
  width = 5
  height = 4
  tmp = data.frame(x = NA, y=NA)
  for (i in seq(6,nrow(amis_mod$mod), 5)){
    tmp = rbind(tmp,c(mean(amis_mod$mod$range[(i-5):i]),log(1/var(amis_mod$mod$logratio[(i-5):i]))))
  }
  pqr_inla_rw2 <- function(data,mod,params,type){
    quants = c(0.1,0.25,0.5,0.75,0.9)
    mu = mod$summary.linear.predictor[,1]
    tau = exp(params[1] + params[2]*data)
    res = data.frame(x = NA, y = NA, quants = NA)
    for (i in seq(length(quants))){
      if (type == "gaussian"){
        tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
      }else if (type == "gamma"){
        tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau
      }
      res = rbind(res,data.frame(x = data, y = tmpy, quants = rep(toString(quants[i]),length(data))))
    }
    return(res[-1,])
  }
  
  pqr_inla_rw2_more <- function(data,mod,params,type){
    quants = c(0.025,0.05,0.15,0.2,0.3,0.35,0.4,0.45,0.55,0.6,0.65,0.7,0.8,0.85,0.95,0.975)
    mu = mod$summary.linear.predictor[,1]
    tau = exp(params[1] + params[2]*data)
    res = data.frame(x = NA, y = NA, quants = NA)
    for (i in seq(length(quants))){
      if (type == "gaussian"){
        tmpy = qnorm(quants[i],mean = mu,sd = 1/sqrt(tau))
      }else if (type == "gamma"){
        tmpy = exp(mu)*qgamma(quants[i], shape = tau, scale = 1)/tau
      }
      res = rbind(res,data.frame(x = data, y = tmpy, quants = rep(toString(quants[i]),length(data))))
    }
    return(res[-1,])
  }
  
  amis_mod$pqr = pqr_inla_rw2(amis_mod$mod$range,res_amis,amis_mod$params,type="gaussian")
  amis_mod$pqr_more = pqr_inla_rw2_more(amis_mod$mod$range,res_amis,amis_mod$params,type="gaussian")
  pqr_text = data.frame(x = c(),y = c(),text = c())
  for (x in unique(amis_mod$pqr$quants)){
    tmpx = amis_mod$pqr$x[amis_mod$pqr$quants==x] + 10
    tmpy = amis_mod$pqr$y[amis_mod$pqr$quants==x]
    tmptext = x
    pqr_text = rbind(pqr_text,data.frame(x=tmpx[length(tmpx)],y=tmpy[length(tmpy)],text = x))
  }
  p6 <- ggplot()+
    geom_line(data = amis_mod$pqr, aes(x=x,y=y,linetype = quants)) +
    geom_line(data = amis_mod$pqr_more, aes(x=x,y=y,color = quants)) +
    geom_point(data = amis_mod$mod, aes(x = range, y = logratio),alpha = 0.4)+
    labs(x="range",y="logratio",color = "", linetype = "") +
    scale_color_manual(values = rep("grey",18)) +
    scale_linetype_manual(values = rep("solid",5)) +
    geom_text(data = pqr_text, aes(x=x,y=y,label= text)) +
    guides(color=F,linetype=F) +
    theme_bw()
  ggsave(filename = "pqr_lidar.pdf", plot = p6, device = NULL, path = "./figures/pqr/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(p6)
}
fig5sup <- figure5sup()
fig5sup
