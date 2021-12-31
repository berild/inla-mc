library(ggplot2)
library(ggpubr)
load("./sims/sem/sem-amis-w-inla.Rdata")
load("./sims/sem/sem-is-w-inla.Rdata")
load("./sims/sem/sem-mcmc-w-inla.Rdata")
load("./sims/sem/sem-mcmc.Rdata")

# Need to add plotting here
# KDE!


gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_temp = gg_color_hue(4)

p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y,linetype="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y,linetype = "MCMC with INLA") ) + 
  geom_line(data = mcmc_mod$margs$`(Intercept)`, aes(x=x,y=y,linetype = "MCMC")) + 
  labs(linetype = "",x=expression(beta[0]),y="") +
  scale_linetype_manual(values = c(1,2,3))+
  theme_bw() +
  theme(legend.position="bottom")
p1


p2 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,linetype="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="MCMC with INLA")) +
  geom_line(data = mcmc_mod$margs$`log(GDPCAP)`, aes(x=x,y=y,linetype = "MCMC")) + 
  labs(linetype = "",x=expression(beta[1]),y="") +
  scale_linetype_manual(values = c(1,2,3))+
  theme_bw() +
  theme(legend.position="bottom")
p2

p3 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,linetype="IS with INLA")) +
  #geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x=x,y=y,color="MCMC with INLA")) +
  geom_line(data = mcmc_mod$margs$tau, aes(x=x,y=y,linetype="MCMC")) +
  labs(linetype = "",x=expression(tau),y="") +
  scale_linetype_manual(values = c(1,2,3))+
  theme_bw() +
  theme(legend.position="bottom")
p3

p4 <- ggplot() +
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, linetype = "AMIS with INLA"),bins = 6,color = "black") +
  geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, linetype = "IS with INLA"),bins = 6,color = "black") +
  #geom_contour(data = mcmc_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) +
  geom_contour(data = mcmc_mod$eta_kern, aes(x = x, y = y, z = z, linetype = "MCMC"),bins = 6,color = "black") +
  labs(linetype = "",x=expression(rho),y=expression(lambda)) +
  scale_linetype_manual(values = c(1,2,3))+
  theme_bw() +
  theme(legend.position="bottom")
p4

ptot <- ggarrange(p1, p2, p3, p4,ncol=2, nrow=2, common.legend = T, legend="bottom",labels = c("a", "b","c","d"),font.label = list(size = 10))
ptot

sp1 <- spplot(turnout["TURNOUT01"],  col.regions= hcl.colors(40,palette = hcl.pals("sequential")[20], alpha = 1, rev = TRUE),colorkey = list(space = "bottom"))
sp1
sp2 <- spplot(turnout["GDPCAP"], col.regions= hcl.colors(40,palette =  hcl.pals("sequential")[32], alpha = 1, rev = T),colorkey = list(space = "bottom"))
sp2
sptot <- ggarrange(sp1, sp2,ncol=2, nrow=1, common.legend = F)
sptot

calc.theta(is_w_inla_mod$theta,is_w_inla_mod$weight,is_w_inla_mod$eta,nrow(is_w_inla_mod$eta),2)


p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$intercept, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data = is_w_inla_mod$margs$intercept, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$intercept, aes(x=x,y=y,color = "MCMC with INLA") ) + 
  #geom_line(data = mcmc_mod$margs$`(Intercept)`, aes(x=x,y=y,color = "MCMC")) + 
  labs(color = "",x=expression(beta[0]),y="") +
  scale_color_manual(values=col_temp[c(1,3)]) + 
  theme_bw() +
  theme(legend.position="bottom")
p1


p2 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data = is_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$GDPCAP, aes(x=x,y=y,color="MCMC with INLA")) +
  #geom_line(data = mcmc_mod$margs$`log(GDPCAP)`, aes(x=x,y=y,color = "MCMC")) + 
  labs(color = "",x=expression(beta[1]),y="") +
  scale_color_manual(values=col_temp[c(1,3)]) +
  theme_bw() +
  theme(legend.position="bottom")
p2

p3 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x=x,y=y,color="AMIS with INLA")) +
  #geom_line(data = is_w_inla_mod$margs$tau, aes(x=x,y=y,color="IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x=x,y=y,color="MCMC with INLA")) +
  #geom_line(data = mcmc_mod$margs$tau, aes(x=x,y=y,color="MCMC")) +
  labs(color = "",x=expression(tau),y="") +
  scale_color_manual(values=col_temp[c(1,3)]) +
  theme_bw() +
  theme(legend.position="bottom")
p3

p4 <- ggplot() +
  geom_contour(data = amis_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "AMIS with INLA"),bins = 6) +
  #geom_contour(data = is_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "IS with INLA"),bins = 6) +
  geom_contour(data = mcmc_w_inla_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC with INLA"),bins = 6) +
  #geom_contour(data = mcmc_mod$eta_kern, aes(x = x, y = y, z = z, color = "MCMC"),bins = 6) +
  labs(color = "",x=expression(rho),y=expression(lambda),linetype="") +
  scale_color_manual(values=col_temp[c(1,3)]) +
  theme_bw() +
  theme(legend.position="bottom")
p4

ptot <- ggarrange(p1, p2, p3, p4,ncol=2, nrow=2, common.legend = T, legend="bottom")
ptot
