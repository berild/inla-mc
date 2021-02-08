library(ggplot2)
library(scales)
library(ggpubr)

## loading simulation results
load(file = "./missing/sims/missing-is-w-inla.Rdata")
load(file = "./missing/sims/missing-amis-w-inla.Rdata")
load(file = "./missing/sims/missing-mcmc-w-inla.Rdata")

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

width = 7
height = 7

p1 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta0, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$beta0, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta0, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(beta[0])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1
ggsave(filename = "missing_beta0.pdf", plot = p1, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
  
p2 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta1, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$beta1, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta1, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(beta[1])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p2
ggsave(filename = "missing_beta1.pdf", plot = p2, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p3 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta2, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$beta2, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta2, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(beta[2])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p3
ggsave(filename = "missing_beta2.pdf", plot = p3, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p4 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$beta3, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$beta3, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$beta3, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(beta[3])) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p4
ggsave(filename = "missing_beta3.pdf", plot = p4, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p5 <- ggplot() + 
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,color="AMIS with INLA")) + 
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x=x,y=y,color="MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(tau)) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5)) + 
  coord_cartesian(xlim=c(0,0.025))
p5
ggsave(filename = "missing_tau.pdf", plot = p5, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)


amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)

p15 <- ggplot() + 
  #geom_hline(yintercept = 10000) + 
  geom_line(data = is_w_inla_mod$ess, aes(x = time, y = ess, color = "IS with INLA"))+
  geom_line(data = amis_w_inla_mod$ess, aes(x = time, y = ess, color = "AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$ess, aes(x = time, y = ess, color = "MCMC with INLA")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() + 
  labs(color = "",x="Runtime (sec)",y="Effective sample size",title="") + 
  #coord_cartesian(xlim=c(min(amis_w_inla_mod$ess$time),max(mcmc_w_inla_mod$ess$time))) + 
  theme_bw() + 
  theme(legend.position="bottom")
p15
ggsave(filename = "missing_ess_plot.pdf", plot = p15, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

ptot1 <- ggarrange(p1, p2, p3, p4, p5, p15, ncol=2, nrow=3, common.legend = TRUE, legend="bottom")
ptot1
ggsave(filename = "missing_tot1.pdf", plot = ptot1, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

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
p6 <-ggplot() + 
  geom_line(data=amis_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[1]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p6
ggsave(filename = "missing_obs1.pdf", plot = p6, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p7 <-ggplot() + 
  geom_line(data=amis_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[2]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[2]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p7
ggsave(filename = "missing_obs2.pdf", plot = p7, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
p8 <-ggplot() + 
  geom_line(data=amis_kerns[[3]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[3]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[3]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[3]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p8
ggsave(filename = "missing_obs3.pdf", plot = p8, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p9 <-ggplot() + 
  geom_line(data=amis_kerns[[4]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[4]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[4]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[4]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p9
ggsave(filename = "missing_obs4.pdf", plot = p9, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p10 <-ggplot() + 
  geom_line(data=amis_kerns[[5]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[5]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[5]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[5]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p10
ggsave(filename = "missing_obs5.pdf", plot = p10, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p11 <-ggplot() + 
  geom_line(data=amis_kerns[[6]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[6]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[6]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[6]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p11
ggsave(filename = "missing_obs6.pdf", plot = p11, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p12 <-ggplot() + 
  geom_line(data=amis_kerns[[7]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[7]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[7]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[7]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p12
ggsave(filename = "missing_obs7.pdf", plot = p12, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p13 <-ggplot() + 
  geom_line(data=amis_kerns[[8]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[8]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[8]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[8]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p13
ggsave(filename = "missing_obs8.pdf", plot = p13, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

p14 <-ggplot() + 
  geom_line(data=amis_kerns[[9]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data=is_kerns[[9]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data=mcmc_kerns[[9]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title=obs.names[9]) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p14
ggsave(filename = "missing_obs9.pdf", plot = p14, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)

ptot2 <- ggarrange(p6, p7, p8, p9, p10, p11, p12, p13, p14, ncol=3, nrow=3, common.legend = TRUE, legend="bottom")
ptot2
ggsave(filename = "missing_tot2.pdf", plot = ptot2, device = NULL, path = "./missing/figures/",
       scale = 1, width = width, height = height, units = "in", dpi=5000)
