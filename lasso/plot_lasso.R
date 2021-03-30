library(ggplot2)
library(scales)
library(ggpubr)

library(ISLR)
library(glmnet)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
col_temp = gg_color_hue(4)

data(Hitters)
Hitters <- na.omit(Hitters)
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


## loading simulation results
load(file = "./sims/lasso/lasso-is.Rdata")
load(file = "./sims/lasso//lasso-amis.Rdata")
load(file = "./sims/lasso/lasso-mcmc.Rdata")
source("./genFuncs.R")

p1 <- ggplot() +
  geom_line(data = amis_w_inla_mod$margs$tau, aes(x = x, y = y, color = "AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$margs$tau, aes(x = x, y = y, color = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$margs$tau, aes(x = x, y = y, color = "MCMC with INLA")) +
  labs(color = "",x="",y="",title=expression(tau)) +
  theme_bw() +
  scale_color_manual(values=col_temp) +
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p1

amis_kerns = lapply(seq(ncol(amis_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = amis_w_inla_mod$eta[,x],
                        weights = amis_w_inla_mod$weight/sum(amis_w_inla_mod$weight),
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})

is_kerns = lapply(seq(ncol(is_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = is_w_inla_mod$eta[,x],
                        bw = 0.08,
                        weights = is_w_inla_mod$weight/sum(is_w_inla_mod$weight),
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})

mcmc_kerns = lapply(seq(ncol(mcmc_w_inla_mod$eta)), function(x){
  as.data.frame(density(x = mcmc_w_inla_mod$eta[,x],
                        from = -0.6, to = 0.6, kernel = "gaussian")[c(1,2)])
})

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)

plot_uni <- function(){
  height2 = 3
  width2 = 5
  height3 = height2
  width3 = width2
  size = 0.8
  tsize = 18
  isize = 18
  asize = 13
  p2 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[2],type = "Lasso Coef."),
               aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = amis_kerns[[1]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_kerns[[1]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_kerns[[1]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="AtBat",y="",linetype = "",title = "a") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.4,0.3)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p3 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[3],type = "Lasso Coef."),
               aes(xintercept = x),color = "salmon3",size = size) +
    geom_line(data = amis_kerns[[2]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_kerns[[2]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_kerns[[2]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="Hits",y="",linetype = "",title = "b") +
    scale_linetype_manual(values = c(1,2,3))+
    theme_bw() +
    coord_cartesian(xlim =c(-0.2,0.6)) +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p4 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[4],type = "Lasso Coef."),
               aes(xintercept = x),color = "salmon3",size = size) +
    geom_line(data = amis_kerns[[3]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_kerns[[3]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_kerns[[3]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="HmRun",y="",linetype = "",title = "c") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.3,0.4)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p5 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[5],type = "Lasso Coef."),
               aes(xintercept = x),color ="salmon3",size = size) +
    geom_line(data = amis_kerns[[4]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_kerns[[4]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_kerns[[4]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="Runs",y="",linetype = "",title = "d") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.25,0.45)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p6 <-ggplot() +
    geom_vline(data = data.frame(x = lasso.coef[6],type = "Lasso Coef."),
               aes(xintercept = x),color = "salmon3",size = size) +
    geom_line(data = amis_kerns[[5]], aes(x=x,y=y,linetype="AMIS with INLA"),size = size) +
    geom_line(data= is_kerns[[5]], aes(x=x,y=y,linetype="IS with INLA"),size = size) +
    geom_line(data= mcmc_kerns[[5]], aes(x=x,y=y,linetype="MCMC with INLA"),size = size) +
    labs(color = "",x="RBI",y="",linetype = "",title = "e") +
    scale_linetype_manual(values = c(1,2,3))+
    coord_cartesian(xlim =c(-0.3,0.6)) +
    theme_bw() +
    theme(legend.position="none",axis.title.y = element_blank(), axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  p7 <- ggplot() +
    geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA"),size = size) +
    geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA"),size = size) +
    geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA"),size = size) +
    scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h","2h"),trans="log",breaks=c(0,60,60*5,60*20,60*60,2*60*60)) +
    labs(linetype = "",color = "",x="Runtime",y="Effective sample size",title = "f") +
    theme_bw() +
    coord_cartesian(xlim = c(10,3*60*70)) +
    scale_linetype_manual(values = c(1,2,3)) +
    theme(legend.position="none",axis.title.y = element_blank(),axis.title.x = element_text(size = tsize),plot.title = element_text(hjust = -0.04, vjust = 0,face="bold",size = isize),axis.text = element_text(size = asize),panel.background = element_blank())
  ggsave(filename = "lasso_uni_atbat.pdf", plot = p2, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_hits.pdf", plot = p3, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_hmrun.pdf", plot = p4, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_runs.pdf", plot = p5, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_rbi.pdf", plot = p6, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  ggsave(filename = "lasso_uni_ess.pdf", plot = p7, device = NULL, path = "./figures/lasso/",
         scale = 1, width = width3, height = height3, units = "in", dpi=5000)
  return(list(p2,p3,p4,p5,p6,p7))
}

plot_uni()

amis_adaptive = lapply(seq(7), function(x){
  tmp = rst(10000,sigma = sqrt(amis_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3,
            mu = amis_w_inla_mod$theta$a.mu[T_s[x],param])
  tmp = sort(tmp)
  y = dst(tmp,mu=amis_w_inla_mod$theta$a.mu[T_s[x]],
          sigma = sqrt(amis_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3)
  y = y/max(y) + x
  data.frame(x=tmp,y=y)
})

truth = lapply(seq(7), function(x){
  tmp = amis_kerns[[param]]$x
  y = amis_kerns[[param]]$y
  y = y/max(y) +x
  data.frame(x=tmp,y=y)
})
is_adaptive = lapply(seq(2), function(x){
  tmp = rst(10000,sigma = sqrt(is_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3,
            mu = is_w_inla_mod$theta$a.mu[T_s[x],param])
  tmp = sort(tmp)
  y = dst(tmp,mu=is_w_inla_mod$theta$a.mu[T_s[x],param],
          sigma = sqrt(is_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3)
  y = y/max(y) + x
  data.frame(x=tmp,y=y)
})

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

res = plot_adaptive(3)
res

cumsum_df <- function(mod){
  weights = mod$weight/sum(mod$weight)
  max(weights)
  n.var <- ncol(mod$eta)
  n.samples = nrow(mod$eta)
  res = list()
  for(i in 1:n.var){
    # Reorder weights
    idx <- order(mod$eta[, i])
    res[[i]] = data.frame(x = 1:n.samples / n.samples, y = cumsum(weights[idx]))
  }
  return(res)
}

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
tmp123 <- cumsum_df2(amis_w_inla_mod)
amiscs <- cumsum_df(amis_w_inla_mod)


cumsumplot <- function(){
  amiscs <- cumsum_df2(amis_w_inla_mod)
  iscs <- cumsum_df2(is_w_inla_mod)
  csamis1 <- ggplot()+
    geom_path(data = amiscs[[4]], aes(x=x,y=y)) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
    labs(x="",y="",title = "AMIS-INLA Runs", linetype ="") +
    theme_bw() +
    theme(plot.title = element_text(size=8,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,11,-8.5), "points"),axis.ticks.length.x = unit(-1, "mm"),
          axis.line = element_line(colour = "gray"),axis.text = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "gray"))

  csamis2 <- ggplot()+
    geom_path(data = amiscs[[5]], aes(x=x,y=y)) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
    labs(x="",y="",title = "AMIS-INLA RBI", linetype ="") +
    theme_bw() +
    theme(plot.title = element_text(size=8,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(-21,5.5,5.5,-8.5), "points"),
          axis.line.y = element_line(colour="gray"),axis.line.x = element_line(colour="black"), axis.text.y = element_blank(),panel.background = element_blank(),panel.border = element_rect(colour = "gray"),plot.background = element_rect(fill = alpha('white',0)))

  csis1 <- ggplot()+
    geom_path(data = iscs[[4]], aes(x=x,y=y)) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
    labs(x="",y="",title = "IS-INLA Runs", linetype ="") +
    theme_bw() +
    theme(plot.title = element_text(size=8,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(5.5,5.5,11,5.5), "points"),axis.ticks.length.x = unit(-1, "mm"),
          axis.line.y = element_line(colour = "black"),axis.line.x = element_line(colour = "gray"), axis.text.x = element_blank(),
          panel.background = element_blank(),panel.border = element_rect(colour="gray"))

  csis2 <- ggplot()+
    geom_path(data = iscs[[5]], aes(x=x,y=y)) +
    geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
    labs(x="",y="",title = "IS-INLA RBI", linetype ="") +
    theme_bw() +
    theme(plot.title = element_text(size=8,vjust=-14,hjust=0.1),axis.title = element_blank(),plot.margin=unit(c(-21,5.5,5.5,5.5), "points"),
          axis.line = element_line(colour="black"), axis.line.x.top = element_line(colour = "white"),panel.background = element_blank(),panel.border = element_rect(colour = "gray"),plot.background = element_rect(fill = alpha('white',0)))

  ptot <- grid.arrange(csis1,csamis1,csis2,csamis2,ncol=2,left=textGrob("Empirical cumulative distribution", rot=90, gp=gpar(fontsize=10),hjust = 0.5),bottom = textGrob("Theoretical cumulative distribution",gp=gpar(fontsize=10),vjust = -0.1))
  return(ptot)
}
cumsumplot()

cumsumplot2 <- function(){
  amiscs = cumsum_df2(amis_w_inla_mod)
  iscs = cumsum_df2(is_w_inla_mod)
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
  return(ptot)
}
cumsumplot2()
