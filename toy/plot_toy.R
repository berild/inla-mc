library(ggplot2)
library(ggpubr)
source("./inlaMC/inlaMC.R")
load("./sims/toy/is_toy.Rdata")
load("./sims/toy/inla_toy.Rdata")
load("./sims/toy/amis_toy.Rdata")
load("./sims/toy/mcmc_toy.Rdata")

# Function to create figure 1
figure1 <- function(param){
  require(LaplacesDemon)
  require(gridExtra)
  require(grid)
  width = 7.5
  height = 3
  T_s = c(1,2,5,10,15,20,28)
  amis_adaptive = lapply(seq(7), function(x){
    tmp = rnorm(500,amis_mod$theta[[1]][T_s[x],1],sqrt(amis_mod$theta[[2]][1,1,T_s[x]]))
    tmp = sort(tmp)
    y = dnorm(tmp,amis_mod$theta[[1]][T_s[x],1],sqrt(amis_mod$theta[[2]][1,1,T_s[x]]))
    y = y/max(y) + x
    data.frame(x=tmp,y=y)
  })
  # Truth are the INLA marginals
  truth = lapply(seq(7), function(x){
    tmp = inla_mod$marginals.fixed$x1[,1]
    y = inla_mod$marginals.fixed$x1[,2]
    y = y/max(y) +x
    data.frame(x=tmp,y=y)
  })
  is_adaptive = lapply(seq(2), function(x){
    tmp = rnorm(500,is_mod$theta[[1]][x,1],sqrt(is_mod$theta[[2]][1,1,x]))
    tmp = sort(tmp)
    y = dnorm(tmp,is_mod$theta[[1]][x,1],sqrt(is_mod$theta[[2]][1,1,x]))
    y = x + y/max(y)
    data.frame(x=tmp,y=y)
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
  ptot <- grid.arrange(pis,pamis,ncol=2,left=textGrob(expression(beta[1]), rot=90, gp=gpar(fontsize=12),hjust = 0.5),bottom = textGrob("Number of adaptations",gp=gpar(fontsize=12),vjust = -0.1))
  ggsave(filename = "toy_adapt.pdf", plot = ptot, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(ptot)
}
fig1 <- figure1(1)
fig1

# Funnction to create Figure 2
figure2 <- function(){
  height = 3
  width = 5
  # Intercept (beta0)
  p1 <- ggplot() +
    geom_line(data = amis_mod$margs$intercept, aes(x=x,y=y, linetype = "A")) +
    geom_line(data = is_mod$margs$intercept, aes(x=x,y=y, linetype = "B")) +
    geom_line(data = mcmc_mod$margs$intercept, aes(x=x,y=y, linetype = "C")) +
    geom_line(data = data.frame(inla_mod$marginals.fixed$`(Intercept)`),aes(x=x,y=y,linetype = "D")) +
    geom_vline(data = data.frame(x = 1), aes(xintercept = 1),color = "salmon3") +
    scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
    labs(color = "", x = expression(beta[0]),y="",linetype = "",title = "a")+
    scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
    theme_bw() +
    coord_cartesian(xlim = c(0,2))+
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
  # Precision (tau)
  p2 <- ggplot() +
    geom_line(data = amis_mod$margs$tau, aes(x=x,y=y,linetype = "A")) +
    geom_line(data = is_mod$margs$tau, aes(x=x,y=y,linetype = "B")) +
    geom_line(data = is_mod$margs$tau, aes(x=x,y=y,linetype = "C")) +
    geom_line(data = data.frame(inla_mod$marginals.hyperpar$`Precision for the Gaussian observations`),aes(x=x,y=y,linetype = "D")) +
    geom_vline(data = data.frame(x = .11), aes(xintercept = 1), color = "salmon3") +
    scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
    scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
    labs(color = "", x = expression(tau),y="",linetype = "",title = "b")+
    theme_bw() +
    coord_cartesian(xlim = c(0.35,2))+
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
  # Fixed effect (beta1)
  p3 <- ggplot() +
    geom_line(data = amis_mod$eta_kern[[1]], aes(x=x,y=y,linetype="A")) +
    geom_line(data = is_mod$eta_kern[[1]], aes(x=x,y=y,linetype="B")) +
    geom_line(data = mcmc_mod$eta_kern[[1]], aes(x=x,y=y,linetype="C")) +
    geom_line(data = data.frame(inla_mod$marginals.fixed$x1),aes(x=x,y=y,linetype = "D")) +
    geom_vline(data = data.frame(x = 1), aes(xintercept = 1), color = "salmon3") +
    labs(color = "", x = expression(beta[1]),y="",linetype = "",title = "c") +
    scale_linetype_manual(values = c(1,2,3,4),labels = c("AMIS-INLA","IS-INLA","MCMC-INLA","INLA")) +
    scale_y_continuous(labels=scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels=scales::number_format(accuracy = 0.1)) +
    coord_cartesian(xlim = c(-0.5,2.5))+
    theme_bw() +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.02,size = 12),axis.title.y = element_blank(), axis.title.x = element_text(size=14))
  # Fixed effect (beta2)
  p4 <- ggplot() +
    geom_line(data = amis_mod$eta_kern[[2]], aes(x=x,y=y,linetype="A")) +
    geom_line(data = is_mod$eta_kern[[2]], aes(x=x,y=y,linetype="B")) +
    geom_line(data = mcmc_mod$eta_kern[[2]], aes(x=x,y=y,linetype="C")) +
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
fig2 <- figure2()
ggarrange(fig2[[1]],fig2[[2]],fig2[[3]],fig2[[4]],ncol=2, nrow=2, common.legend = T,legend="bottom")

# Calculating the running effective sample sizes
amis_mod$ess = running.ESS(amis_mod$eta,amis_mod$times,ws = amis_mod$weight,norm = T)
is_mod$ess = running.ESS(is_mod$eta,is_mod$times,ws = is_mod$weight,norm = T)
mcmc_mod$ess = running.ESS(mcmc_mod$eta,mcmc_mod$times) # THIS TAKES TIME!

# kernel density estimation of the joint distribution of beta1 and beta2
eta_joint_kern_amis = kde2d.weighted(x = amis_mod$eta[,1], y = amis_mod$eta[,2], w = amis_mod$weight/(sum(amis_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_is = kde2d.weighted(x = is_mod$eta[,1], y = is_mod$eta[,2], w = is_mod$weight/(sum(is_mod$weight)), n = 100, lims = c(0,2,-2,0))
eta_joint_kern_mcmc = kde2d.weighted(x = mcmc_mod$eta[,1], y = mcmc_mod$eta[,2], w = rep(1/length(mcmc_mod$eta[,1]),length(mcmc_mod$eta[,1])), n = 100, lims = c(0,2,-2,0))
amis_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_amis$x, y=eta_joint_kern_amis$y), z=as.vector(eta_joint_kern_amis$z))
is_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_is$x, y=eta_joint_kern_is$y), z=as.vector(eta_joint_kern_is$z))
mcmc_mod$eta_joint_kern = data.frame(expand.grid(x=eta_joint_kern_mcmc$x, y=eta_joint_kern_mcmc$y), z=as.vector(eta_joint_kern_mcmc$z))

# Function to create Figure 3
figure3 <- function(){
  nbin=12
  size = 0.8
  height = 3
  width = 5
  cont1 <- ggplot() +
    geom_contour(data = amis_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="a") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))
  
  cont2 <- ggplot() +
    geom_contour(data = amis_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="b") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))
  
  cont3 <- ggplot() +
    geom_contour(data = amis_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "AMIS with INLA"),bins = nbin,color="black") +
    geom_contour(data = is_mod$eta_joint_kern, aes(x = x+100, y = y+100, z = z, linetype = "IS with INLA"),bins = nbin,color="black") +
    geom_contour(data = mcmc_mod$eta_joint_kern, aes(x = x, y = y, z = z, linetype = "MCMC with INLA"),bins = nbin,color="black") +
    labs(linetype="",color = "",x=expression(beta[1]),y=expression(beta[2]),linetype="",title="c") +
    geom_point(data = data.frame(x = 1,y = -1), aes(x=x,y=y), shape = 4, stroke = 2,size = 3,color = "salmon3") +
    theme_bw() +
    scale_linetype_manual(values = c(1,2,3)) +
    coord_cartesian(xlim = c(0,2),ylim=c(-2,0)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))
  
  essp <- ggplot() +
    geom_line(data = amis_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA")) +
    geom_line(data = is_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA")) +
    geom_line(data = mcmc_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA")) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title="d") +
    theme_bw() +
    coord_cartesian(xlim = c(10,60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",plot.title = element_text(face = "bold",hjust = -0.04,size = 12),axis.title = element_text(size=14))
  ggsave(filename = "toy_joint_amis.pdf", plot = cont1, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "toy_joint_is.pdf", plot = cont2, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "toy_joint_mcmc.pdf", plot = cont3, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  ggsave(filename = "toy_joint_ess.pdf", plot = essp, device = NULL, path = "./figures/toy/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(list(cont1,cont2,cont3,essp))
}
fig3<- figure3()
ggarrange(fig3[[1]],fig3[[2]],fig3[[3]],fig3[[4]],ncol=2, nrow=2, common.legend = T,legend="bottom")
