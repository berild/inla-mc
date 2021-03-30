library(ggplot2)
library(scales)
library(ggpubr)
library(mice) # data

## loading simulation results
load(file = "./sims/missing/missing-is.Rdata")
load(file = "./sims/missing//missing-amis.Rdata")
load(file = "./sims/missing/missing-mcmc.Rdata")

data(nhanes2)

d.mis <- nhanes2
idx.mis <- which(is.na(d.mis$bmi))
n.mis <- length(idx.mis)

df = list(d.mis = d.mis, idx.mis = idx.mis)

marginal.stat <- function(method){
  require(INLA)
  stats = lapply(method, function(y){
    m1 = inla.emarginal(function(x) x,y)
    m2 = inla.emarginal(function(x) x^2,y)
    c(m1,sqrt(m2-m1^2))
  })
  return(stats)
}

amis_w_inla_mod$stats = marginal.stat(amis_w_inla_mod$margs)
is_w_inla_mod$stats = marginal.stat(is_w_inla_mod$margs)
mcmc_w_inla_mod$stats = marginal.stat(mcmc_w_inla_mod$margs)

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)

plot_uni <- function(){
  height2 = 3
  width2 = 5
  height3 = height2
  width3 = width2
  size = 0.8
  tsize = 19
  isize = 18
  asize = 13
  p1 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$beta0, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$margs$beta0, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$margs$beta0, aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    scale_linetype_manual(values = c(1,2,3)) +
    labs(linetype = "",x=expression(beta[0]),y="",title = "a") +
    coord_cartesian(xlim = c(-200,300))+
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = isize),
          panel.background = element_blank(),axis.text = element_text(size = asize))
  p2 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$beta1, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$margs$beta1, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$margs$beta1, aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(linetype = "",y="",x=expression(beta[1]),title = "b") +
    coord_cartesian(xlim = c(-5,15))+
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = isize),
          panel.background = element_blank(),axis.text = element_text(size = asize))
  p3 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$beta2, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$margs$beta2, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$margs$beta2, aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(linetype = "",y="",x=expression(beta[2]),title="c") +
    coord_cartesian(xlim = c(-50,100))+
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = isize),
          panel.background = element_blank(),axis.text = element_text(size = asize))
  p4 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$beta3, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$margs$beta3, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$margs$beta3, aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(linetype = "",y="",x=expression(beta[3]),title = "d") +
    coord_cartesian(xlim = c(-50,150))+
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = isize),
          panel.background = element_blank(),axis.text = element_text(size = asize))
  p5 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",y="",x=expression(tau),title = "e") +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = isize),
          panel.background = element_blank(),axis.text = element_text(size = asize)) +
    coord_cartesian(xlim=c(0,0.025))
  essp <- ggplot() +
    geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA"),size = size) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h","2h"),trans="log",breaks=c(0,60,60*5,60*20,60*60,2*60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title = "f") +
    theme_bw() +
    coord_cartesian(xlim = c(10,2*60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.02, vjust = 0,face="bold",size = isize),
          panel.background = element_blank(),axis.text = element_text(size = asize))
  ggsave(filename = "missing_uni_beta0.pdf", plot = p1, device = NULL, path = "./figures/missing/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "missing_uni_beta1.pdf", plot = p2, device = NULL, path = "./figures/missing/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "missing_uni_beta2.pdf", plot = p3, device = NULL, path = "./figures/missing/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "missing_uni_beta3.pdf", plot = p4, device = NULL, path = "./figures/missing/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "missing_uni_tau.pdf", plot = p5, device = NULL, path = "./figures/missing/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "missing_uni_ess.pdf", plot = essp, device = NULL, path = "./figures/missing/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  return(list(p1,p2,p3,p4,p5,essp))
}

plot_uni()

amis_kerns = lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
is_kerns = lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
mcmc_kerns = lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        kernel = "gaussian")[c(1,2)])
})
plot_eta <- function(){
  obs.names = sprintf("Observation %d",df$idx.mis)
  height2 = 3
  width2 = 5
  height3 = height2
  width3 = width2
  size = 0.8
  isize = 15
  asize = 13
  res = list()
  for (i in seq(length(obs.names))){
    ptmp <- ggplot() +
      geom_line(data=amis_kerns[[i]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
      geom_line(data=is_kerns[[i]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
      geom_line(data=mcmc_kerns[[i]], aes(x=x,y=y,linetype="MCMC with INLA"),size =size) +
      labs(linetype = "",x="",y="",title=obs.names[i]) +
      scale_linetype_manual(values = c(1,2,3))+
      theme_bw() +
      theme(legend.position="none",axis.title = element_blank(),plot.title = element_text(hjust = 0.5, vjust = 0,face="bold",size = isize),
           panel.background = element_blank(),axis.text = element_text(size = asize)) +
      coord_cartesian(xlim = c(0,60))
    ggsave(filename = sprintf("missing_obs%d.pdf",i), plot = ptmp, device = NULL, path = "./figures/missing/",
           scale = 1, width = width2, height = height2, units = "in", dpi=5000)
    res[[i]] = ptmp
  }
  return(res)
}
plot_eta()
