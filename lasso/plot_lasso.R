library(ggplot2)
library(scales)
library(ggpubr)

library(ISLR)
library(glmnet)
source("./inlaMC/inlaMC.R")
load("./sims/lasso/mcmc_lasso.Rdata")
load("./sims/lasso/is_lasso.Rdata")
load("./sims/lasso/amis_lasso.Rdata")


data(Hitters)
Hitters <- na.omit(Hitters)

# GETTING ML estimates
#Create variables for lasso
x <- model.matrix(Salary ~ ., Hitters)[, -1]
x <- x[, 1:5] #Just for testing
x <- scale(x)
y <- Hitters$Salary
y <- scale(y)
df <- list(y = y, x = x)

#Indices for train/test model
set.seed(1)
train <- sample(1:nrow(x), nrow(x)/2)
test <- (-train)


#Grid for lambda parameter in lasso
grid <- 10^seq(10, -2, length = 100)

#Fit model to complete dataset
out <- glmnet(x, y, alpha = 1, lambda = grid,intercept=F)
lasso.coef <- predict(out, type = "coefficients", s = 0.073,intercept = F)
lasso.coef


amis_mod$eta_kern = kde_mc(amis_mod$eta, amis_mod$weight)
is_mod$eta_kern = kde_mc(is_mod$eta, is_mod$weight)
mcmc_mod$eta_kern = kde_mc(mcmc_mod$eta,rep(1,nrow(mcmc_mod$eta)))


# need to change the binwidth for the kde in IS-INLA to not only getting point masses
is_mod$eta_kern = lapply(seq(ncol(is_mod$eta)), function(x){
  as.data.frame(density(x = is_mod$eta[,x],
                        bw = 0.08,
                        weights = is_mod$weight/sum(is_mod$weight),
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})


# Calculating the running effective sanple size
amis_mod$ess = running.ESS(amis_mod$eta, amis_mod$times,ws =  amis_mod$weight)
is_mod$ess = running.ESS(is_mod$eta, is_mod$times,ws =  is_mod$weight)
mcmc_mod$ess = running.ESS(mcmc_mod$eta, mcmc_mod$times)

# Function to create figure 4 in the paper
figure4 <- function(){
  height2 = 3
  width2 = 5
  height3 = height2
  width3 = width2
  size = 0.8
  tsize = 18
  isize = 18
  asize = 13
  p1 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[2],type = "Lasso Coef."),
               aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = amis_mod$eta_kern[[1]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_mod$eta_kern[[1]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_mod$eta_kern[[1]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="AtBat",y="",linetype = "",title = "a") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p2 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[3],type = "Lasso Coef."),
               aes(xintercept = x),color = "salmon3",size = size) +
    geom_line(data = amis_mod$eta_kern[[2]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_mod$eta_kern[[2]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_mod$eta_kern[[2]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="Hits",y="",linetype = "",title = "b") +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    coord_cartesian(xlim =c(-0.2,0.6)) +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  
  p3 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[4],type = "Lasso Coef."),
               aes(xintercept = x),color = "salmon3",size = size) +
    geom_line(data = amis_mod$eta_kern[[3]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_mod$eta_kern[[3]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_mod$eta_kern[[3]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="HmRun",y="",linetype = "",title = "c") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.3,0.4)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p4 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[5],type = "Lasso Coef."),
               aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = amis_mod$eta_kern[[4]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_mod$eta_kern[[4]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_mod$eta_kern[[4]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="Runs",y="",linetype = "",title = "d") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.25,0.45)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p5 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[6],type = "Lasso Coef."),
               aes(xintercept = x),color = "salmon3",size = size) +
    geom_line(data = amis_mod$eta_kern[[5]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_mod$eta_kern[[5]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_mod$eta_kern[[5]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="RBI",y="",linetype = "",title = "e") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.3,0.6)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p6 <- ggplot() +
    geom_line(data = amis_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA"),size = size) +
    geom_line(data = is_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA"),size = size) +
    geom_line(data = mcmc_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA"),size = size) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h","2h"),trans="log",breaks=c(0,60,60*5,60*20,60*60,2*60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title = "f") +
    theme_bw() +
    coord_cartesian(xlim = c(10,3*60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  ggsave(filename = "lasso_uni_atbat.pdf", plot = p1, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_hits.pdf", plot = p2, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_hmrun.pdf", plot = p3, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_runs.pdf", plot = p4, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_rbi.pdf", plot = p5, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_ess.pdf", plot = p6, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  return(list(p1,p2,p3,p4,p5,p6))
}

fig4 <- figure4()
ggarrange(fig4[[1]],fig4[[2]],fig4[[3]],fig4[[4]],fig4[[5]],fig4[[6]],ncol=3, nrow=2, common.legend = T,legend="bottom")


figure5 <- function(param){
  require(LaplacesDemon)
  require(gridExtra)
  require(grid)
  width = 7.5
  height = 3
  T_s = c(1,2,5,10,15,20,28)
  amis_adaptive = lapply(seq(7), function(x){
    tmpx = seq(-1,1,length.out = 1000)
    tmpy = dst(tmpx,mu=amis_mod$theta[[1]][T_s[x]],
            sigma = sqrt(amis_mod$theta[[2]][param,param,T_s[x]]), nu=3)
    tmpy = tmpy/max(tmpy) + x
    data.frame(x = tmpx,y = tmpy)
  })

  truth = lapply(seq(7), function(x){
    tmp = amis_mod$eta_kern[[param]]$x
    y = amis_mod$eta_kern[[param]]$y
    y = y/max(y) +x
    data.frame(x=tmp,y=y)
  })
  is_adaptive = lapply(seq(2), function(x){
    tmpx = seq(-1,1,length.out = 1000)
    tmpy = dst(tmpx,mu=is_mod$theta[[1]][T_s[x],param],
              sigma = sqrt(is_mod$theta[[2]][param,param,T_s[x]]), nu=3)
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
  ggsave(filename = "lasso_adapt.pdf", plot = ptot, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(ptot)
}
figure5(3)

# Function to calculate probability plots of the weights
cumsum_df2 <- function(mod){
  weights = mod$weight/sum(mod$weight)
  max(weights)
  n.var <- ncol(mod$eta)
  n.samples = nrow(mod$eta)
  res = list()
  for(i in 1:n.var){
    # Reorder weights
    idx <- order(mod$eta[, i])
    tmp = data.frame(x = 1:n.samples / n.samples, y = cumsum(weights[idx]))
    tmp_s = spline(tmp$x,tmp$y,n = 500)
    res[[i]] = as.data.frame(tmp_s)
  }
  return(res)
}


# Function to create figure 6 in the paper
figure6 <- function(){
  height = 4.5
  width = 7.5
  amiscs = cumsum_df2(amis_mod)
  iscs = cumsum_df2(is_mod)
  csamis1 <- ggplot()+
    geom_line(data = amiscs[[4]], aes(x=x,y=y),linetype = 1) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),color = "salmon3") +
    labs(x="",y="",title = "Runs", linetype ="",color="") +
    theme_bw() +
    theme(plot.title = element_text(size=10,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,11,-8.5), "points"),axis.ticks.length.x = unit(-1, "mm"),
          axis.line = element_line(colour = "gray"),axis.text = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "gray"))

  csamis2 <- ggplot()+
    geom_line(data = amiscs[[5]], aes(x=x,y=y),linetype = 1) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),color = "salmon3") +
    labs(x="",y="",title = "RBI", linetype ="",color="") +
    theme_bw() +
    theme(plot.title = element_text(size=10,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(-23,5.5,5.5,-8.5), "points"),
          axis.line.y = element_line(colour="gray"),axis.line.x = element_line(colour="black"), axis.text.y = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "gray"),plot.background = element_rect(fill = alpha('white',0)))

  csis1 <- ggplot()+
    geom_line(data = iscs[[4]], aes(x=x,y=y),linetype = 2) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),color = "salmon3") +
    labs(x="",y="",title = "Runs", linetype ="",color="") +
    theme_bw() +
    theme(plot.title = element_text(size=10,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,11,5.5), "points"),axis.ticks.length.x = unit(-1, "mm"),
          axis.line.y = element_line(colour = "black"),axis.line.x = element_line(colour = "gray"), axis.text.x = element_blank(),
          panel.background = element_blank(),panel.border = element_rect(colour="gray"))

  csis2 <- ggplot()+
    geom_line(data = iscs[[5]], aes(x=x,y=y),linetype = 2) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),color = "salmon3") +
    labs(x="",y="",title = "RBI", linetype ="",color="") +
    theme_bw() +
    theme(plot.title = element_text(size=10,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(-23,5.5,5.5,5.5), "points"),
          axis.line = element_line(colour="black"), axis.line.x.top = element_line(colour = "white"),panel.background = element_blank(),panel.border = element_rect(colour = "gray"),plot.background = element_rect(fill = alpha('white',0)))
  ptot <- grid.arrange(csis1,csamis1,csis2,csamis2,ncol=2,left=textGrob("Empirical cumulative distribution", rot=90, gp=gpar(fontsize=12),hjust = 0.5),bottom = textGrob("Theoretical cumulative distribution",gp=gpar(fontsize=12),vjust = -0.1))
  ggsave(filename = "lasso_cumsum2.pdf", plot = ptot, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width, height = height, units = "in", dpi=5000)
  return(ptot)
}
figure6()
