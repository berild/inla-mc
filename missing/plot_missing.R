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

p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$beta0, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$beta0, aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta0, aes(x=x,y=y,linetype="MCMC with INLA")) +
  scale_linetype_manual(values = c(1,2,3)) +
  labs(linetype = "",x=expression(beta[0]),y="") +
  coord_cartesian(xlim = c(-200,300))+
  theme_bw() +
  theme(legend.position="none",axis.title = element_text(size = 12))
p1

p2 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$beta1, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$beta1, aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta1, aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(linetype = "",y="",x=expression(beta[1])) +
  coord_cartesian(xlim = c(-5,15))+
  theme_bw() +
  theme(legend.position="none",axis.title = element_text(size = 12))
p2

p3 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$beta2, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$beta2, aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta2, aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(linetype = "",y="",x=expression(beta[2])) +
  coord_cartesian(xlim = c(-50,100))+
  theme_bw() +
  theme(legend.position="none",axis.title = element_text(size = 12))
p3
p4 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$beta3, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$beta3, aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta3, aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(linetype = "",y="",x=expression(beta[3])) +
  coord_cartesian(xlim = c(-50,150))+
  theme_bw() +
  theme(legend.position="none",axis.title = element_text(size = 12))
p4

p5 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(color = "",y="",x=expression(tau)) +
  theme_bw() +
  theme(legend.position="none",axis.title = element_text(size = 12))+
  coord_cartesian(xlim=c(0,0.025))
p5


amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)

essp <- ggplot() +
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA")) +
  scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) +
  labs(linetype = "",color = "",x="Runtime",y="Effective sample size") +
  theme_bw() +
  coord_cartesian(xlim = c(10,60*70)) +
  scale_linetype_manual(values = c(1,2,3)) +
  theme(legend.position="none",axis.title = element_text(size = 12))
essp

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

obs.names = sprintf("Observation %d",df$idx.mis)
for (i in seq(length(obs.names))){
  ptmp <- ggplot() +
    geom_line(data=amis_kerns[[i]], aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data=is_kerns[[i]], aes(x=x,y=y,linetype="IS with INLA")) +
    geom_line(data=mcmc_kerns[[i]], aes(x=x,y=y,linetype="MCMC with INLA")) +
    labs(linetype = "",x="",y="",title=obs.names[i]) +
    theme_bw() +
    theme(legend.position="none",plot.title = element_text(size = 12,vjust=-8,hjust=0.1)) +
    coord_cartesian(xlim = c(0,60))
  ggsave(filename = sprintf("missing_obs%d.pdf",i), plot = ptmp, device = NULL, path = "./figures/missing/",
         scale = 1, width = 5, height = 3, units = "in", dpi=5000)
}

