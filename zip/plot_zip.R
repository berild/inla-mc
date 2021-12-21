library(ggplot2)
library(pscl)

source("genFuncs.R")

load("./zip/is_zip.Rdata")
load("./zip/amis_zip.Rdata")
height2 = 3
width2 = 5
height3 = height2
width3 = width2
size = 0.8
tsize = 18
isize = 18
asize = 13

zinb <- read.csv("./zip/fish.csv")
zinb <- within(zinb, {
  nofish <- factor(nofish)
  livebait <- factor(livebait)
  camper <- factor(camper)
})

summary(zinb)


# Model:
#  Poisson: count ~ child + camper 
# Binomial: ~ 1 + persons

#summary(m1 <- zeroinfl(count ~ child + camper | persons, data = zinb))


# Reduce data for testing
#d <- zinb[1:30, ]

# Full data
d <- zinb

# ZIP
ml.res <- summary(zeroinfl(count ~ child + camper | persons, data = d))
ml.res

amis_mod$ess = running.ESS(amis_mod$eta, amis_mod$times,ws =  amis_mod$weight)
is_mod$ess = running.ESS(is_mod$eta, is_mod$times,ws =  is_mod$weight)
plot_uni <- function(){
  p1 <-ggplot() +
    geom_vline(data = data.frame(x = ml.res$coefficients$count[1,1],type = "ML estimate"),aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = is_mod$margs$intercept, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = amis_mod$margs$intercept, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    labs(color = "",x="Intercept",y="",linetype = "",title = "a") +
    scale_linetype_manual(values = c(1,2,3))+
    #coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),
          plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),
          axis.text = element_text(size = asize),panel.background = element_blank())
  
  p2 <-ggplot() +
    geom_vline(data = data.frame(x = ml.res$coefficients$count[2,1],type = "ML estimate"),aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = is_mod$margs$child, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = amis_mod$margs$child, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    labs(color = "",x="Child",y="",linetype = "",title = "b") +
    scale_linetype_manual(values = c(1,2,3))+
    #coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),
          plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),
          axis.text = element_text(size = asize),panel.background = element_blank())
  
  
  p3 <-ggplot() +
    geom_vline(data = data.frame(x = ml.res$coefficients$count[3,1],type = "ML estimate"),aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = is_mod$margs$camper1, aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = amis_mod$margs$camper1, aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    labs(color = "",x="Camper1",y="",linetype = "",title = "c") +
    scale_linetype_manual(values = c(1,2,3))+
    #coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),
          plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),
          axis.text = element_text(size = asize),panel.background = element_blank())
  
  
  p4 <-ggplot() +
    geom_vline(data = data.frame(x = ml.res$coefficients$zero[1,1],type = "ML estimate"),aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = is_mod$eta_kern[[1]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = amis_mod$eta_kern[[1]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    labs(color = "",x=expression(alpha[1]),y="",linetype = "",title = "d") +
    scale_linetype_manual(values = c(1,2,3))+
    #coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),
          plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),
          axis.text = element_text(size = asize),panel.background = element_blank())
  
  p5 <-ggplot() +
    geom_vline(data = data.frame(x = ml.res$coefficients$zero[2,1],type = "ML estimate"),aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = is_mod$eta_kern[[2]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data = amis_mod$eta_kern[[2]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    labs(color = "",x=expression(alpha[2]),y="",linetype = "",title = "e") +
    scale_linetype_manual(values = c(1,2,3))+
    #coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),
          plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),
          axis.text = element_text(size = asize),panel.background = element_blank())
  
  p6 <- ggplot() +
    geom_line(data = amis_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA"),size = size) +
    geom_line(data = is_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA"),size = size) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h","2h"),trans="log",breaks=c(0,60,60*5,60*20,60*60,2*60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title = "f") +
    theme_bw() +
    coord_cartesian(xlim = c(10,3*60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  
  ggsave(filename = "zip_uni_intercept.pdf", plot = p1, device = NULL, path = "./zip/figures/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "zip_uni_child.pdf", plot = p2, device = NULL, path = "./zip/figures/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "zip_uni_camper1.pdf", plot = p3, device = NULL, path = "./zip/figures/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "zip_uni_alpha1.pdf", plot = p4, device = NULL, path = "./zip/figures/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "zip_uni_alpha2.pdf", plot = p5, device = NULL, path = "./zip/figures/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "zip_uni_ess.pdf", plot = p6, device = NULL, path = "./zip/figures/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  return(list(p1,p2,p3,p4,p5))
}

