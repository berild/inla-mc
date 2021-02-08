library(ggplot2)
library(ggpubr)
load("./toy/toy-is.Rdata")
load("./toy/toy-inla.Rdata")
load("./toy/toy-amis.Rdata")
load("./toy/toy-mcmc.Rdata")
source("./toy/genFuncs.R")


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_temp = gg_color_hue(4)

p1 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y, color = "MCMC with INLA")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$`(Intercept)`),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1, linetype = "Truth")) + 
  scale_linetype_manual(values= "dotdash") + 
  scale_color_manual(values = col_temp) + 
  labs(color = "", x = expression(beta[0]),y="",linetype = "")+
  theme_bw() + 
  coord_cartesian(xlim = c(-0.5,2.5))+
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1

p2 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y, color = "AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y, color = "IS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y, color = "MCMC with INLA")) +
  geom_line(data = data.frame(inla_mod$marginals.hyperpar$`Precision for the Gaussian observations`),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1, linetype = "Truth")) + 
  labs(color = "", x = expression(tau),y="",linetype = "")+
  scale_linetype_manual(values= "dotdash") + 
  theme_bw() + 
  scale_color_manual(values = col_temp) + 
  coord_cartesian(xlim = c(0.3,2.5))+
  theme(legend.position="bottom")
p2

amis_w_inla_mod$eta_uni_kerns = lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
is_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        kernel = "gaussian")[c(1,2)])
})
mcmc_w_inla_mod$eta_uni_kerns= lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        weights = rep(1/length(mcmc_w_inla_mod$eta[,x]),length(mcmc_w_inla_mod$eta[,x])),
                        kernel = "gaussian")[c(1,2)])
})

eta_joint_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1/length(mcmc_w_inla_mod$eta[,1]),length(mcmc_w_inla_mod$eta[,1])), n = 100, lims = c(0,2,-2,0))
amis_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))
is_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_is$x, y=eta_joint_kern_is$y), z=as.vector(eta_joint_kern_is$z))
mcmc_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_mcmc$x, y=eta_joint_kern_mcmc$y), z=as.vector(eta_joint_kern_mcmc$z))

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight, norm = F)
save(amis_w_inla_mod,file="./sims/toy/toy-amis-w-inla.Rdata")
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight, norm = F)
save(is_w_inla_mod,file="./sims/toy/toy-is-w-inla.Rdata")
mcmc_w_inla_mod$ess = running.ESS(mcmc_w_inla_mod$eta, mcmc_w_inla_mod$times, norm = F)
save(mcmc_w_inla_mod,file="./sims/toy/toy-mcmc-w-inla.Rdata")


p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) + 
  geom_line(data = data.frame(inla_mod$marginals.fixed$x1),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1, linetype = "Truth")) + 
  scale_linetype_manual(values= "dotdash") + 
  labs(color = "", x = expression(beta[1]),y="",linetype = "") + 
  scale_color_manual(values =col_temp) + 
  coord_cartesian(xlim = c(-0.5,2.5))+
  theme_bw() + 
  theme(legend.position="bottom")
p3
p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,color="MCMC with INLA")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$x2),aes(x=x,y=y,color = "INLA")) + 
  geom_vline(data = data.frame(x = 1), aes(xintercept = -1, linetype = "Truth")) + 
  scale_linetype_manual(values= "dotdash") + 
  labs(color = "", x = expression(beta[2]),y="",linetype = "") + 
  scale_color_manual(values = col_temp) + 
  coord_cartesian(xlim = c(-2.5,0.5))+
  theme_bw() + 
  theme(legend.position="bottom")
p4

ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = T,legend="bottom",labels = c("a", "b","c","d"),font.label = list(size = 10))


cont1 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, color = "IS with INLA"),bins = 6) +
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, color = "MCMC with INLA"),bins = 6) +
  labs(color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="") +
  theme_bw() +
  scale_color_manual(values = col_temp[-2]) + 
  coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) + 
  theme(legend.position="bottom")
cont1      

cont2 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, color = "MCMC with INLA"),bins = 6) +
  labs(color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="") +
  theme_bw() +
  scale_color_manual(values = col_temp[-2]) + 
  coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) + 
  theme(legend.position="bottom")
cont2

cont3 <- ggplot() + 
  geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, color = "AMIS with INLA"),bins = 6) +
  geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, color = "IS with INLA"),bins = 6) +
  geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) +
  labs(color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="") +
  theme_bw() +
  scale_color_manual(values = col_temp[-2]) + 
  coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) + 
  theme(legend.position="bottom")
cont3

essp <- ggplot() + 
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,color = "AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,color = "MCMC with INLA")) +
  scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) + 
  labs(color = "",x="Runtime",y="ESS") +
  theme_bw() + 
  coord_cartesian(xlim = c(10,60*70)) + 
  scale_color_manual(values = col_temp[-2]) + 
  theme(legend.position="bottom")
essp

ggarrange(cont1,cont2,cont3,essp,ncol=2, nrow=2, common.legend = T,legend="bottom",labels = c("a", "b","c","d"),font.label = list(size = 10))


p1mcmc <- ggplot() + 
  geom_path(data = data.frame(x = seq(length(mcmc_w_inla_mod$full_eta[1:500,1])), y = mcmc_w_inla_mod$full_eta[1:500,1]),aes(x=x,y=y),color="darkgrey") + 
  geom_path(data = data.frame(x = seq(501,500+length(mcmc_w_inla_mod$eta[,1])), y = mcmc_w_inla_mod$eta[,1]),aes(x=x,y=y)) + 
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  labs(x = "N",y = expression(beta[1]),color="",fill="") +
  coord_cartesian(ylim = c(-0.5,2.5)) + 
  theme_bw()
p1mcmc

p2mcmc <- ggplot() + 
  geom_path(data = data.frame(x = seq(length(mcmc_w_inla_mod$full_eta[1:500,2])), y = mcmc_w_inla_mod$full_eta[1:500,2]),aes(x=x,y=y),color = "darkgrey") + 
  geom_path(data = data.frame(x = seq(501,500 + length(mcmc_w_inla_mod$eta[,2])), y = mcmc_w_inla_mod$eta[,2]),aes(x=x,y=y)) +
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  labs(color = "",fill = "") + 
  labs(x = "N",y = expression(beta[2]),color="",fill="") +
  coord_cartesian(ylim = c(-2.5,0.5)) + 
  theme_bw()
p2mcmc

ggarrange(p1mcmc,p2mcmc,ncol=2, nrow=1, common.legend = T,legend="bottom")

T_s = c(1,2,5,10,15,20,28)
amis_adaptive = lapply(seq(7), function(x){
  tmp = rnorm(500,amis_w_inla_mod$theta$a.mu[T_s[x],1],sqrt(amis_w_inla_mod$theta$a.cov[1,1,T_s[x]]))
  tmp = sort(tmp)
  y = dnorm(tmp,amis_w_inla_mod$theta$a.mu[T_s[x],1],sqrt(amis_w_inla_mod$theta$a.cov[1,1,T_s[x]]))
  y = y/max(y) + x
  data.frame(x=tmp,y=y)
})
amis_adaptive2 = lapply(seq(7), function(x){
  tmp = inla_mod$marginals.fixed$x1[,1]
  y = inla_mod$marginals.fixed$x1[,2]
  y = y/max(y) +x
  data.frame(x=tmp,y=y)
})

p1amis <- ggplot() + 
  geom_polygon(data = amis_adaptive2[[1]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[2]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[3]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[4]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[5]],aes(x=y,y=x,fill = "target")) +
  geom_polygon(data = amis_adaptive2[[6]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive2[[7]],aes(x=y,y=x,fill = "target")) + 
  geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) + 
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  coord_cartesian(ylim = c(-6,6))+
  scale_x_continuous(label = T_s - 1, breaks = seq(7)) + 
  labs(x = "t",y = expression(beta[1]),color="",fill="",title="AMIS with INLA") +
  theme_bw() + 
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
p1amis

amis_adaptive21 = lapply(seq(7), function(x){
  tmp = rnorm(500,amis_w_inla_mod$theta$a.mu[T_s[x],2],sqrt(amis_w_inla_mod$theta$a.cov[2,2,T_s[x]]))
  tmp = sort(tmp)
  y = dnorm(tmp,amis_w_inla_mod$theta$a.mu[T_s[x],2],sqrt(amis_w_inla_mod$theta$a.cov[2,2,T_s[x]]))
  y = y/max(y) + x
  data.frame(x=tmp,y=y)
})
amis_adaptive22 = lapply(seq(7), function(x){
  tmp = inla_mod$marginals.fixed$x2[,1]
  y = inla_mod$marginals.fixed$x2[,2]
  y = y/max(y) +x
  data.frame(x=tmp,y=y)
})

p2amis <-ggplot() + 
  geom_polygon(data = amis_adaptive22[[1]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[2]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[3]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[4]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[5]],aes(x=y,y=x,fill = "target")) +
  geom_polygon(data = amis_adaptive22[[6]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = amis_adaptive22[[7]],aes(x=y,y=x,fill = "target")) + 
  geom_path(data = amis_adaptive21[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[2]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[3]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[4]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[5]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[6]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = amis_adaptive21[[7]],aes(x=y,y=x, color = "proposal")) +
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  coord_cartesian(ylim = c(-6,6))+
  scale_x_continuous(label = T_s-1, breaks = seq(7)) + 
  labs(x = "t",y = expression(beta[2]),color="",fill="",title="AMIS with INLA") +
  theme_bw() + 
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
p2amis

is_adaptive = lapply(seq(2), function(x){
  tmp = rnorm(500,is_w_inla_mod$theta$a.mu[x,1],sqrt(is_w_inla_mod$theta$a.cov[1,1,x]))
  tmp = sort(tmp)
  y = dnorm(tmp,is_w_inla_mod$theta$a.mu[x,1],sqrt(is_w_inla_mod$theta$a.cov[1,1,x]))
  y = x + y/max(y)
  data.frame(x=tmp,y=y)
})
is_adaptive2 = lapply(seq(2), function(x){
  tmp = inla_mod$marginals.fixed$x1[,1]
  y = inla_mod$marginals.fixed$x1[,2]
  y = x + y/max(y)
  data.frame(x=tmp,y=y)
})
p1is <-ggplot() + 
  geom_polygon(data = is_adaptive2[[1]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = is_adaptive2[[2]],aes(x=y,y=x,fill = "target")) + 
  geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal")) + 
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  scale_x_continuous(label = seq(2)-1, breaks = seq(2)) + 
  coord_cartesian(ylim = c(-6,6))+
  labs(color="",fill="",x="t",y = expression(beta[1]),title="IS with INLA") +
  theme_bw()+
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
p1is

is_adaptive21 = lapply(seq(2), function(x){
  tmp = rnorm(500,is_w_inla_mod$theta$a.mu[x,2],sqrt(is_w_inla_mod$theta$a.cov[2,2,x]))
  tmp = sort(tmp)
  y = dnorm(tmp,is_w_inla_mod$theta$a.mu[x,2],sqrt(is_w_inla_mod$theta$a.cov[2,2,x]))
  y = x + y/max(y)
  data.frame(x=tmp,y=y)
})
is_adaptive22 = lapply(seq(2), function(x){
  tmp = inla_mod$marginals.fixed$x2[,1]
  y = inla_mod$marginals.fixed$x2[,2]
  y = x + y/max(y)
  data.frame(x=tmp,y=y)
})

p2is <-ggplot() + 
  geom_polygon(data = is_adaptive22[[1]],aes(x=y,y=x,fill = "target")) + 
  geom_polygon(data = is_adaptive22[[2]],aes(x=y,y=x,fill = "target")) + 
  geom_path(data = is_adaptive21[[1]],aes(x=y,y=x, color = "proposal")) + 
  geom_path(data = is_adaptive21[[2]],aes(x=y,y=x, color = "proposal")) + 
  scale_color_manual(values = col_temp[1]) + 
  scale_fill_manual(values = col_temp[3]) + 
  scale_x_continuous(label = seq(2)-1, breaks = seq(2)) + 
  labs(color="",fill="",x="t",y = expression(beta[2]),title="IS with INLA") +
  coord_cartesian(ylim = c(-6,6))+
  theme_bw()+ 
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
p2is

ggarrange(p1amis,p1is,ncol=2, nrow=1, common.legend = T,legend="bottom",labels = c("a", "b"),font.label = list(size = 10))
