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



p2 <-ggplot() +
  geom_vline(data = data.frame(x = lasso.coef[2],type = "Lasso Coef."),
             aes(xintercept = x),color ="salmon3") +
  geom_line(data = amis_kerns[[1]], aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data= is_kerns[[1]], aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data= mcmc_kerns[[1]], aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(color = "",x="AtBat",y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3))+
  coord_cartesian(xlim =c(-0.4,0.3)) +
  theme_bw() +
  theme(legend.position="none")
p2

p3 <-ggplot() +
  geom_vline(data = data.frame(x = lasso.coef[3],type = "Lasso Coef."),
             aes(xintercept = x),color = "salmon3") +
  geom_line(data = amis_kerns[[2]], aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data= is_kerns[[2]], aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data= mcmc_kerns[[2]], aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(color = "",x="Hits",y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3))+
  theme_bw() +
  coord_cartesian(xlim =c(-0.2,0.6)) +
  theme(legend.position="none")
p3

p4 <-ggplot() +
  geom_vline(data = data.frame(x = lasso.coef[4],type = "Lasso Coef."),
             aes(xintercept = x),color = "salmon3") +
  geom_line(data = amis_kerns[[3]], aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data= is_kerns[[3]], aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data= mcmc_kerns[[3]], aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(color = "",x="HmRun",y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3))+
  coord_cartesian(xlim =c(-0.3,0.4)) +
  theme_bw() +
  theme(legend.position="none")
p4

p5 <-ggplot() +
  geom_vline(data = data.frame(x = lasso.coef[5],type = "Lasso Coef."),
             aes(xintercept = x),color ="salmon3") +
  geom_line(data = amis_kerns[[4]], aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data= is_kerns[[4]], aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data= mcmc_kerns[[4]], aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(color = "",x="Runs",y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3))+
  coord_cartesian(xlim =c(-0.25,0.45)) +
  theme_bw() +
  theme(legend.position="none")
p5

p6 <-ggplot() +
  geom_vline(data = data.frame(x = lasso.coef[6],type = "Lasso Coef."),
             aes(xintercept = x),color = "salmon3") +
  geom_line(data = amis_kerns[[5]], aes(x=x,y=y,linetype="AMIS with INLA")) +
  geom_line(data= is_kerns[[5]], aes(x=x,y=y,linetype="IS with INLA")) +
  geom_line(data= mcmc_kerns[[5]], aes(x=x,y=y,linetype="MCMC with INLA")) +
  labs(color = "",x="RBI",y="",linetype = "") +
  scale_linetype_manual(values = c(1,2,3))+
  coord_cartesian(xlim =c(-0.3,0.6)) +
  theme_bw() +
  theme(legend.position="none")
p6

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)


p7 <- ggplot() +
  geom_line(data = amis_w_inla_mod$ess,aes(x=time,y=ess,linetype = "AMIS with INLA")) +
  geom_line(data = is_w_inla_mod$ess,aes(x=time,y=ess,linetype = "IS with INLA")) +
  geom_line(data = mcmc_w_inla_mod$ess,aes(x=time,y=ess,linetype = "MCMC with INLA")) +
  scale_x_continuous(labels = c("0 sec", "1 min", "5 min", "20 min", "1 h"),trans="log",breaks=c(0,60,60*5,60*20,60*60)) +
  labs(linetype = "",color = "",x="Runtime",y="Effective sample size") +
  theme_bw() +
  coord_cartesian(xlim = c(10,60*70)) +
  scale_linetype_manual(values = c(1,2,3)) +
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))
p7

ptot <- ggarrange(p2, p3, p4, p5, p6, p7, ncol=2, nrow=3, common.legend = TRUE, legend="bottom",labels = c("a","b","c","d","e","f"),font.label = list(size = 10))
ptot



plot_adaptive <- function(param){
  require(LaplacesDemon)
  T_s = c(1,2,5,10,15,20,28)
  amis_adaptive = lapply(seq(7), function(x){
    tmp = rst(1000,sigma = sqrt(amis_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3,
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
    tmp = rst(1000,sigma = sqrt(is_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3,
              mu = is_w_inla_mod$theta$a.mu[T_s[x],param])
    tmp = sort(tmp)
    y = dst(tmp,mu=is_w_inla_mod$theta$a.mu[T_s[x],param],
            sigma = sqrt(is_w_inla_mod$theta$a.cov[param,param,T_s[x]]), nu=3)
    y = y/max(y) + x
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
    coord_cartesian(ylim = c(-0.5,0.5))+
    scale_x_continuous(label = T_s - 1, breaks = seq(7)) +
    labs(x = "Number of adaptations",y = "HmRun", color="",fill="",title="AMIS-INLA") +
    theme_bw() +
    theme(plot.title = element_text(size=8,vjust=-8,hjust=0.1),panel.grid.minor = element_blank(),legend.position = "none",panel.grid.major.y = element_blank())
  pis <-ggplot() +
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target"),alpha = 0.7) +
    geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal"),linetype=2) +
    geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal"),linetype=2) +
    scale_color_manual(values = "black") +
    scale_fill_manual(values = "salmon3") +
    scale_x_continuous(label = seq(2)-1, breaks = seq(2)) +
    coord_cartesian(xlim = c(0,3.5),ylim = c(-0.5,0.5))+
    labs(color="",fill="",x="Number of adaptations",y = "HmRun",title="IS-INLA") +
    theme_bw()+
    theme(plot.title = element_text(size=8,vjust=-8,hjust=0.1),panel.grid.minor = element_blank(),legend.position = "none",panel.grid.major.y = element_blank())
  return(list(pis,pamis))
}

res = plot_adaptive(3)
res[[1]]
res[[2]]

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

amiscs <- cumsum_df(amis_w_inla_mod)

csamis1 <- ggplot()+
  geom_path(data = amiscs[[4]], aes(x=x,y=y)) +
  geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
  labs(x="Theoretical cumulative distribution",y="Empirical cumulative distribution",title = "AMIS-INLA Runs", linetype ="") +
  theme_bw() +
  theme(plot.title = element_text(size=12,vjust=-8,hjust=0.1))
csamis1

csamis2 <- ggplot()+
  geom_path(data = amiscs[[5]], aes(x=x,y=y)) +
  geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
  labs(x="Theoretical cumulative distribution",y="Empirical cumulative distribution",title = "AMIS-INLA RBI", linetype ="") +
  theme_bw() +
  theme(plot.title = element_text(size=12,vjust=-8,hjust=0.1))
csamis2

iscs <- cumsum_df(is_w_inla_mod)

csis1 <- ggplot()+
  geom_path(data = iscs[[4]], aes(x=x,y=y)) +
  geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
  labs(x="Theoretical cumulative distribution",y="Empirical cumulative distribution",title = "IS-INLA Runs", linetype ="") +
  theme_bw() +
  theme(plot.title = element_text(size=12,vjust=-8,hjust=0.1))
csis1

csis2 <- ggplot()+
  geom_path(data = iscs[[5]], aes(x=x,y=y)) +
  geom_line(data= data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y),linetype = 2) +
  labs(x="Theoretical cumulative distribution",y="Empirical cumulative distribution",title = "IS-INLA RBI", linetype ="") +
  theme_bw() +
  theme(plot.title = element_text(size=12,vjust=-8,hjust=0.1))
csis2
