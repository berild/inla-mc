library(ggplot2)
library(ggpubr)
source("./inlaMC/inlaMC.R")
load("./sims/toy/is_toy.Rdata")
load("./sims/toy/inla_toy.Rdata")
load("./sims/toy/amis_toy.Rdata")
load("./sims/toy/mcmc_toy.Rdata")

p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y, linetype = "A")) +
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y, linetype = "B")) +
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y, linetype = "C")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$`(Intercept)`),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1),color = "salmon3") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  labs(color = "", x = expression(beta[0]),y="",linetype = "")+
  theme_bw() +
  coord_cartesian(xlim = c(0,2))+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))
p1

p2 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,linetype = "A")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype = "B")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype = "C")) +
  geom_line(data = data.frame(inla_mod$marginals.hyperpar$`Precision for the Gaussian observations`),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = .11), aes(xintercept = 1), color = "salmon3") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  labs(color = "", x = expression(tau),y="",linetype = "")+
  theme_bw() +
  coord_cartesian(xlim = c(0.35,2))+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))
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

p3 <- ggplot() +
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,linetype="A")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,linetype="B")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,linetype="C")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$x1),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1), color = "salmon3") +
  labs(color = "", x = expression(beta[1]),y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  coord_cartesian(xlim = c(-0.5,2.5))+
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))
p3
p4 <- ggplot() +
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,linetype="A")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,linetype="B")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,linetype="C")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$x2),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = 1), aes(xintercept = -1), color = "salmon3") +
  labs(color = "", x = expression(beta[2]),y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  coord_cartesian(xlim = c(-2.5,0.5))+
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))
p4

ggarrange(p1,p2,p3,p4,ncol=2, nrow=2, common.legend = T,legend="bottom",labels = c("a", "b","c","d"),font.label = list(size = 10))

joint_plot <- function(){
  nbin=12
  size = 0.8
  height = 3
  width = 5
  cont1 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="a") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 14),axis.title = element_text(size=14))

  cont2 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="b") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 14),axis.title = element_text(size=14))

  cont3 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="c") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 14),axis.title = element_text(size=14))

  essp <- ggplot() +
    geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA")) +
    geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA")) +
    geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA")) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title="d") +
    theme_bw() +
    coord_cartesian(xlim = c(10,60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 14),axis.title = element_text(size=14))
  ggsave(filename = sprintf("toy_joint_amis.pdf",i), plot = cont1, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = sprintf("toy_joint_is.pdf",i), plot = cont2, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = sprintf("toy_joint_mcmc.pdf",i), plot = cont3, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = sprintf("toy_joint_ess.pdf",i), plot = essp, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(list(cont1,cont2,cont3,essp))
}
joint_plot()
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

plot_adaptive <- function(param){
  require(LaplacesDemon)
  require(gridExtra)
  require(grid)
  T_s = c(1,2,5,10,15,20,28)
  amis_adaptive = lapply(seq(7), function(x){
    tmpx = seq(-1,1,length.out = 1000)
    tmpy = dst(tmpx,mu=amis_w_inla_mod$theta$a.mu[T_s[x]],
               sigma = sqrt(amis_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3)
    tmpy = tmpy/max(tmpy) + x
    data.frame(x = tmpx,y = tmpy)
  })

  truth = lapply(seq(7), function(x){
    tmp = amis_kerns[[param]]$x
    y = amis_kerns[[param]]$y
    y = y/max(y) +x
    data.frame(x=tmp,y=y)
  })
  is_adaptive = lapply(seq(2), function(x){
    tmpx = seq(-1,1,length.out = 1000)
    tmpy = dst(tmpx,mu=is_w_inla_mod$theta$a.mu[T_s[x],param],
               sigma = sqrt(is_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3)
    tmpy = tmpy/max(tmpy) + x
    data.frame(x = tmpx,y = tmpy)
  })
  pamis <- ggplot() +
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[3]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[4]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[5]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[6]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[7]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "salmon3") +
    coord_cartesian(ylim = c(-0.5,0.5))+
    scale_x_continuous(label = T_s - 1, breaks = seq(7)) +
    labs(x = "",y = "", color="",fill="",title="AMIS-INLA") +
    theme_bw() +
    theme(plot.title = element_text(size=10,vjust=-10,hjust=0.06),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,5.5,-8.5), "points"),axis.text.y = element_blank(),
          panel.background = element_blank(),panel.border = element_rect(colour = "gray"),panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),axis.text = element_text(size = 10),legend.position = "none")
  pis <-ggplot() +
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
    geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal"),linetype=2) +
    geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal"),linetype=2) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "salmon3") +
    scale_x_continuous(label = seq(2)-1, breaks = seq(2)) +
    coord_cartesian(xlim = c(0,3.5),ylim = c(-0.5,0.5))+
    labs(color="",fill="",x="",y = "",title="IS-INLA") +
    theme_bw()+
    theme(plot.title = element_text(size=10,vjust=-10,hjust=0.06),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,5.5,5.5), "points"),legend.position = "none",
          panel.background = element_blank(),panel.border = element_rect(colour = "gray"),panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),axis.text = element_text(size = 10))
  ptot <- grid.arrange(pis,pamis,ncol=2,left=textGrob("HmRun", rot=90, gp=gpar(fontsize=12),hjust = 0.5),bottom = textGrob("Number of adaptations",gp=gpar(fontsize=12),vjust = -0.1))
  return(ptot)
}
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
  geom_polygon(data = amis_adaptive2[[1]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[2]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[3]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[4]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[5]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[6]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[7]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "salmon3") +
  coord_cartesian(ylim = c(-6,6))+
  scale_x_continuous(label = T_s - 1, breaks = seq(7)) +
  labs(x = "Number of adaptations",y = expression(beta[1]),color="",fill="",title="AMIS-INLA") +
  theme_bw() +
  theme(plot.title = element_text(size=8,vjust=-8,hjust=0.1),panel.grid.minor = element_blank(),legend.position = "none",panel.grid.major.y = element_blank())
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
  labs(x = "Adaptations",y = expression(beta[2]),color="",fill="",title="AMIS with INLA") +
  theme_bw() +
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.2),panel.grid.minor = element_blank())
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
  geom_polygon(data = is_adaptive2[[1]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
  geom_polygon(data = is_adaptive2[[2]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
  geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal"),linetype=2) +
  geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal"),linetype=2) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "salmon3") +
  scale_x_continuous(label = seq(2)-1, breaks = seq(2)) +
  coord_cartesian(ylim = c(-6,6),xlim=c(0,3.5))+
  labs(color="",fill="",x="Number of adaptations",y = expression(beta[1]),title="IS-INLA") +
  theme_bw()+
  theme(plot.title = element_text(size=8,vjust=-8,hjust=0.1),panel.grid.minor = element_blank(),legend.position = "none",panel.grid.major.y = element_blank())
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
  coord_cartesian(ylim = c(-6,6),xlim=c(-2,4))+
  theme_bw()+
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
p2is

ggarrange(p1amis,p1is,ncol=2, nrow=1, common.legend = T,legend="bottom",labels = c("a", "b"),font.label = list(size = 10))

library(ggplot2)
library(ggpubr)
load("./sims/toy/toy-is.Rdata")
load("./sims/toy/toy-inla.Rdata")
load("./sims/toy/toy-amis.Rdata")
load("./sims/toy/toy-mcmc.Rdata")
source("./genFuncs.R")


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

plot_uni <- function(){
  height = 3
  width = 5
p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y, linetype = "A")) +
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y, linetype = "B")) +
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y, linetype = "C")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$`(Intercept)`),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1),color = "salmon3") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  labs(color = "", x = expression(beta[0]),y="",linetype = "",title = "a")+
  scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
  theme_bw() +
  coord_cartesian(xlim = c(0,2))+
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p2 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,linetype = "A")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype = "B")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype = "C")) +
  geom_line(data = data.frame(inla_mod$marginals.hyperpar$`Precision for the Gaussian observations`),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = .11), aes(xintercept = 1), color = "salmon3") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
  labs(color = "", x = expression(tau),y="",linetype = "",title = "b")+
  theme_bw() +
  coord_cartesian(xlim = c(0.35,2))+
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p3 <- ggplot() +
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,linetype="A")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,linetype="B")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[1]], aes(x=x,y=y,linetype="C")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$x1),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = 1), aes(xintercept = 1), color = "salmon3") +
  labs(color = "", x = expression(beta[1]),y="",linetype = "",title = "c") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
  coord_cartesian(xlim = c(-0.5,2.5))+
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
p4 <- ggplot() +
  geom_line(data = amis_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,linetype="A")) +
  geom_line(data = is_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,linetype="B")) +
  geom_line(data = mcmc_w_inla_mod$eta_uni_kerns[[2]], aes(x=x,y=y,linetype="C")) +
  geom_line(data = data.frame(inla_mod$marginals.fixed$x2),aes(x=x,y=y,linetype = "D")) +
  geom_vline(data = data.frame(x = 1), aes(xintercept = -1), color = "salmon3") +
  labs(color = "", x = expression(beta[2]),y="",linetype = "",title = "d") +
  scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
  scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
  scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
  coord_cartesian(xlim = c(-2.5,0.5))+
  theme_bw() +
  theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
ggsave(filename = "toy_uni_intercept.pdf", plot = p1, device = NULL, path = "./figures/toy/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
ggsave(filename = "toy_uni_eta1.pdf", plot = p3, device = NULL, path = "./figures/toy/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
ggsave(filename = "toy_uni_eta2.pdf", plot = p4, device = NULL, path = "./figures/toy/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
ggsave(filename = "toy_uni_tau.pdf", plot = p2, device = NULL, path = "./figures/toy/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
return(list(p1,p2,p3,p4))
}

plot_uni()

eta_joint_kern_amis = kde2d.weighted(x = amis_w_inla_mod$eta[,1], y = amis_w_inla_mod$eta[,2], w = amis_w_inla_mod$weight/(sum(amis_w_inla_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_is = kde2d.weighted(x = is_w_inla_mod$eta[,1], y = is_w_inla_mod$eta[,2], w = is_w_inla_mod$weight/(sum(is_w_inla_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_mcmc = kde2d.weighted(x = mcmc_w_inla_mod$eta[,1], y = mcmc_w_inla_mod$eta[,2], w = rep(1/length(mcmc_w_inla_mod$eta[,1]),length(mcmc_w_inla_mod$eta[,1])), n = 100, lims = c(0,2,-2,0))
amis_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))
is_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_is$x, y=eta_joint_kern_is$y), z=as.vector(eta_joint_kern_is$z))
mcmc_w_inla_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_mcmc$x, y=eta_joint_kern_mcmc$y), z=as.vector(eta_joint_kern_mcmc$z))

joint_plot <- function(){
  nbin=12
  size = 0.8
  height = 3
  width = 5
  cont1 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="a") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))

  cont2 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="b") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))

  cont3 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_w_inla_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_w_inla_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="c") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))

  essp <- ggplot() +
    geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA")) +
    geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA")) +
    geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA")) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title="d") +
    theme_bw() +
    coord_cartesian(xlim = c(10,60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))
  ggsave(filename = sprintf("toy_joint_amis.pdf",i), plot = cont1, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = sprintf("toy_joint_is.pdf",i), plot = cont2, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = sprintf("toy_joint_mcmc.pdf",i), plot = cont3, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = sprintf("toy_joint_ess.pdf",i), plot = essp, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(list(cont1,cont2,cont3,essp))
}
joint_plot()

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

plot_adaptive <- function(param){
  require(LaplacesDemon)
  require(gridExtra)
  require(grid)
  T_s = c(1,2,5,10,15,20,28)
  amis_adaptive = lapply(seq(7), function(x){
    tmpx = seq(-6,6,length.out = 1000)
    tmpy = dnorm(tmpx,amis_w_inla_mod$theta$a.mu[T_s[x]],sqrt(amis_w_inla_mod$theta$a.cov[param,param,T_s[x]]))
    tmpy = tmpy/max(tmpy) + x
    data.frame(x = tmpx,y = tmpy)
  })
  truth = lapply(seq(7), function(x){
    tmp = inla_mod$marginals.fixed$x1[,1]
    y = inla_mod$marginals.fixed$x1[,2]
    y = y/max(y) +x
    data.frame(x=tmp,y=y)
  })
  is_adaptive = lapply(seq(2), function(x){
    tmpx = seq(-6,6,length.out = 1000)
    tmpy = dnorm(tmpx,is_w_inla_mod$theta$a.mu[T_s[x],param],sqrt(is_w_inla_mod$theta$a.cov[param,param,T_s[x]]))
    tmpy = tmpy/max(tmpy) + x
    data.frame(x = tmpx,y = tmpy)
  })
  pamis <- ggplot() +
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[3]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[4]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[5]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[6]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_polygon(data = truth[[7]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
    geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) +
    geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "salmon3") +
    coord_cartesian(ylim = c(-6,6))+
    scale_x_continuous(label = T_s - 1, breaks = seq(7)) +
    labs(x = "",y = "", color="",fill="",title="AMIS-INLA") +
    theme_bw() +
    theme(plot.title = element_text(size=10,vjust=-10,hjust=0.06),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,5.5,-8.5), "points"),axis.text.y = element_blank(),
          panel.background = element_blank(),panel.border = element_rect(colour = "gray"),panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),axis.text = element_text(size = 10),legend.position = "none")
  pis <-ggplot() +
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
    geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal"),linetype=2) +
    geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal"),linetype=2) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "salmon3") +
    scale_x_continuous(label = seq(2)-1, breaks = seq(2)) +
    coord_cartesian(xlim = c(0,3.5),ylim = c(-6,6))+
    labs(color="",fill="",x="",y = "",title="IS-INLA") +
    theme_bw()+
    theme(plot.title = element_text(size=10,vjust=-10,hjust=0.06),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,5.5,5.5), "points"),legend.position = "none",
          panel.background = element_blank(),panel.border = element_rect(colour = "gray"),panel.grid.minor = element_blank(),panel.grid.major.y = element_blank(),axis.text = element_text(size = 10))
  ptot <- grid.arrange(pis,pamis,ncol=2,left=textGrob("HmRun", rot=90, gp=gpar(fontsize=12),hjust = 0.5),bottom = textGrob("Number of adaptations",gp=gpar(fontsize=12),vjust = -0.1))
  return(ptot)
}
plot_adaptive(1)
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
  geom_polygon(data = amis_adaptive2[[1]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[2]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[3]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[4]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[5]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[6]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_polygon(data = amis_adaptive2[[7]],aes(x=y,y=x,fill = "target"),alpha=0.7) +
  geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) +
  geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "salmon3") +
  coord_cartesian(ylim = c(-6,6))+
  scale_x_continuous(label = T_s - 1, breaks = seq(7)) +
  labs(x = "Number of adaptations",y = expression(beta[1]),color="",fill="",title="AMIS-INLA") +
  theme_bw() +
  theme(plot.title = element_text(size=8,vjust=-8,hjust=0.1),panel.grid.minor = element_blank(),legend.position = "none",panel.grid.major.y = element_blank())
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
  labs(x = "Adaptations",y = expression(beta[2]),color="",fill="",title="AMIS with INLA") +
  theme_bw() +
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.2),panel.grid.minor = element_blank())
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
  geom_polygon(data = is_adaptive2[[1]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
  geom_polygon(data = is_adaptive2[[2]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
  geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal"),linetype=2) +
  geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal"),linetype=2) +
  scale_color_manual(values = "black") +
  scale_fill_manual(values = "salmon3") +
  scale_x_continuous(label = seq(2)-1, breaks = seq(2)) +
  coord_cartesian(ylim = c(-6,6),xlim=c(0,3.5))+
  labs(color="",fill="",x="Number of adaptations",y = expression(beta[1]),title="IS-INLA") +
  theme_bw()+
  theme(plot.title = element_text(size=8,vjust=-8,hjust=0.1),panel.grid.minor = element_blank(),legend.position = "none",panel.grid.major.y = element_blank())
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
  coord_cartesian(ylim = c(-6,6),xlim=c(-2,4))+
  theme_bw()+
  theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
p2is

ggarrange(p1amis,p1is,ncol=2, nrow=1, common.legend = T,legend="bottom",labels = c("a", "b"),font.label = list(size = 10))
