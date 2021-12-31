library(ggplot2)
library(ggpubr)
source("./inlaMC/inlaMC.R")
load("./sims/sem/amis_sem.Rdata")
load("./sims/sem/is_sem.Rdata")
load("./sims/sem/mcmc_sem.Rdata")


# Calculating the joint distribution of rho and lambda (kernel density estimation)
eta_kern_amis = kde2d.weighted(x = amis_mod$eta[,1], y = amis_mod$eta[,2], w = amis_mod$weight/sum(amis_mod$weight), n = 100, lims = c(-1,1,-1,1))
amis_mod$eta_joint_kern = data.frame(expand.grid(x=eta_kern_amis$x, y=eta_kern_amis$y), z=as.vector(eta_kern_amis$z))

eta_kern_is = kde2d.weighted(x = is_mod$eta[,1], y = is_mod$eta[,2], w =is_mod$weight/sum(is_mod$weight), n = 100, lims = c(-1,1,-1,1))
is_mod$eta_joint_kern = data.frame(expand.grid(x=eta_kern_is$x, y=eta_kern_is$y), z=as.vector(eta_kern_is$z))

eta_kern_mc = kde2d.weighted(x = mcmc_mod$samples$rho, y = mcmc_mod$samples$lambda, w = rep(1,nrow(mcmc_mod$samples))/nrow(mcmc_mod$samples), n = 100, lims = c(-1,1,-1,1))
mcmc_mod$eta_joint_kern = data.frame(expand.grid(x=eta_kern_mc$x, y=eta_kern_mc$y), z=as.vector(eta_kern_mc$z))

# Function creating figure 7 in the paper
figure7 <- function(){
  p1 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y,linetype="IS with INLA")) +
    geom_line(data = mcmc_mod$margs$`(Intercept)`, aes(x=x,y=y,linetype = "MCMC")) + 
    labs(linetype = "",x=expression(beta[0]),y="") +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    theme(legend.position="bottom")
  p2 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data = is_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,linetype="IS with INLA")) +
    geom_line(data = mcmc_mod$margs$`log(GDPCAP)`, aes(x=x,y=y,linetype = "MCMC")) + 
    labs(linetype = "",x=expression(beta[1]),y="") +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    theme(legend.position="bottom")
  p3 <- ggplot() +
    geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="AMIS with INLA")) +
    geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="IS with INLA")) +
    geom_line(data = mcmc_mod$margs$tau, aes(x=x,y=y,linetype="MCMC")) +
    labs(linetype = "",x=expression(tau),y="") +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    theme(legend.position="bottom")
  p4 <- ggplot() +
    geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, linetype = "AMIS with INLA"),bins = 6,color = "black") +
    geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, linetype = "IS with INLA"),bins = 6,color = "black") +
    geom_contour(data = mcmc_mod$eta_kern, aes(x = x, y = y, z = z, linetype = "MCMC"),bins = 6,color = "black") +
    labs(linetype = "",x=expression(rho),y=expression(lambda)) +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    theme(legend.position="bottom")
  return(list(p1,p2,p3,p4))
}

fig7 <- figure7()
ptot <- ggarrange(fig7[[1]],fig7[[2]],fig7[[3]],fig7[[4]],ncol=2, nrow=2, labels = c("a","b","c","d"),common.legend = T,legend="bottom")
ptot
ggsave(filename = "sem_tot.pdf", plot = ptot, device = NULL, path = "./figures/sem/",
       scale = 1, width = 4.5, height = 4.5, units = "in", dpi=5000)
