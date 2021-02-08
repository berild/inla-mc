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
load(file = "./sims/lasso/lasso-is-w-inla.Rdata")
load(file = "./sims/lasso//lasso-amis-w-inla.Rdata")
load(file = "./sims/lasso/lasso-mcmc-w-inla.Rdata")
source("./lasso/lasso_general_function.R")

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
             aes(xintercept = x,linetype = type)) + 
  geom_line(data = amis_kerns[[1]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[1]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[1]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="AtBat",linetype = "") + 
  scale_linetype_manual(values = "dashed" )+
  scale_color_manual(values=col_temp) + 
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p2

p3 <-ggplot() + 
  geom_vline(data = data.frame(x = lasso.coef[3],type = "Lasso Coef."), 
             aes(xintercept = x,linetype = type)) + 
  geom_line(data = amis_kerns[[2]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[2]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[2]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="Hits",linetype = "") + 
  scale_linetype_manual(values = "dashed" )+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p3

p4 <-ggplot() + 
  geom_vline(data = data.frame(x = lasso.coef[4],type = "Lasso Coef."), 
             aes(xintercept = x,linetype = type)) + 
  geom_line(data = amis_kerns[[3]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[3]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[3]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="HmRun",linetype = "") + 
  scale_linetype_manual(values = "dashed" )+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p4

p5 <-ggplot() + 
  geom_vline(data = data.frame(x = lasso.coef[5],type = "Lasso Coef."), 
             aes(xintercept = x,linetype = type)) +  
  geom_line(data = amis_kerns[[4]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[4]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[4]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="Runs",linetype = "") + 
  scale_linetype_manual(values = "dashed" )+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p5

p6 <-ggplot() + 
  geom_vline(data = data.frame(x = lasso.coef[6],type = "Lasso Coef."), 
             aes(xintercept = x,linetype = type)) +  
  geom_line(data = amis_kerns[[5]], aes(x=x,y=y,color="AMIS with INLA")) +
  geom_line(data= is_kerns[[5]], aes(x=x,y=y,color="IS with INLA")) + 
  geom_line(data= mcmc_kerns[[5]], aes(x=x,y=y,color="MCMC with INLA")) + 
  labs(color = "",x="",y="",title="RBI",linetype = "") + 
  scale_linetype_manual(values = "dashed" )+
  theme_bw() + 
  theme(legend.position="bottom",plot.title = element_text(hjust = 0.5))
p6

amis_w_inla_mod$ess = running.ESS(amis_w_inla_mod$eta, amis_w_inla_mod$times,ws =  amis_w_inla_mod$weight)
is_w_inla_mod$ess = running.ESS(is_w_inla_mod$eta, is_w_inla_mod$times,ws =  is_w_inla_mod$weight)

p7 <- ggplot() + 
  #geom_hline(yintercept = 10000) + 
  geom_line(data = is_w_inla_mod$ess, aes(x = time, y = ess, color = "IS with INLA"))+
  geom_line(data = amis_w_inla_mod$ess, aes(x = time, y = ess, color = "AMIS with INLA")) + 
  geom_line(data = mcmc_w_inla_mod$ess, aes(x = time, y = ess, color = "MCMC with INLA")) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  annotation_logticks() + 
  labs(color = "",x="Runtime (sec)",y="Effective sample size",title="") + 
  coord_cartesian(xlim=c(min(amis_w_inla_mod$ess$time),max(mcmc_w_inla_mod$ess$time))) + 
  theme_bw() + 
  theme(legend.position="bottom")
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
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target")) + 
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target")) + 
    geom_polygon(data = truth[[3]],aes(x=y,y=x,fill = "target")) + 
    geom_polygon(data = truth[[4]],aes(x=y,y=x,fill = "target")) + 
    geom_polygon(data = truth[[5]],aes(x=y,y=x,fill = "target")) +
    geom_polygon(data = truth[[6]],aes(x=y,y=x,fill = "target")) + 
    geom_polygon(data = truth[[7]],aes(x=y,y=x,fill = "target")) + 
    geom_path(data = amis_adaptive[[1]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = amis_adaptive[[2]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = amis_adaptive[[3]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = amis_adaptive[[4]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = amis_adaptive[[5]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = amis_adaptive[[6]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = amis_adaptive[[7]],aes(x=y,y=x, color = "proposal")) + 
    scale_color_manual(values = col_temp[1]) + 
    scale_fill_manual(values = col_temp[3]) + 
    coord_cartesian(ylim = c(-0.5,0.5))+
    scale_x_continuous(label = T_s - 1, breaks = seq(7)) + 
    labs(x = "t",y = "HmRun",color="",fill="",title="AMIS with INLA") +
    theme_bw() + 
    theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
  pis <-ggplot() + 
    geom_polygon(data = truth[[1]],aes(x=y,y=x,fill = "target")) + 
    geom_polygon(data = truth[[2]],aes(x=y,y=x,fill = "target")) + 
    geom_path(data = is_adaptive[[1]],aes(x=y,y=x, color = "proposal")) + 
    geom_path(data = is_adaptive[[2]],aes(x=y,y=x, color = "proposal")) + 
    scale_color_manual(values = col_temp[1]) + 
    scale_fill_manual(values = col_temp[3]) + 
    scale_x_continuous(label = seq(2)-1, breaks = seq(2)) + 
    coord_cartesian(ylim = c(-0.5,0.5))+
    labs(color="",fill="",x="t",y = "HmRun",title="IS with INLA") +
    theme_bw()+
    theme(plot.title = element_text(size=7,vjust=-8,hjust=0.01))
  ptot <- ggarrange(pamis,pis,ncol=2, nrow=1, common.legend = T,legend="bottom",labels = c("a", "b"),font.label = list(size = 10))
  return(ptot)
}

plot_adaptive(3)
