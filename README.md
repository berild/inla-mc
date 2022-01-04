# Importance Sampling with the Integrated Nested Laplace Approximation

This repository contains the code used in the Importance Sampling with the Integrated Nested Laplace Approximation paper. The implementation of the Importance Sampling with the Integrated Nested Laplace Approximation (IS-INLA) can be found in <a href="https://github.com/berild/inla-mc/blob/master/inlaMC/inlaIS.R">inlaIS.R</a>, the Adaptive Multiple Importance Sampling with the Integrated Nested Laplace Approximation (AMIS-INLA) in <a href="https://github.com/berild/inla-mc/blob/master/inlaMC/inlaAMIS.R">inlaAMIS.R</a>, and the Markov Chain Monte Carlo with the Integrated Nested Laplace Approximation (MCMC-INLA) in <a href="https://github.com/berild/inla-mc/blob/master/inlaMC/inlaMH.R">inlaMH.R</a>. In addition, general functions used in all algorithms are collected in <a href="https://github.com/berild/inla-mc/blob/master/inlaMC/genFuncs.R">genFuncs.R</a>. 

In this repository, we have included the code for the respective examples presented in the paper and supplementary materials. These include

#### Paper
* Bivariate Linear Gaussian model (<a href="https://github.com/berild/inla-mc/blob/master/toy/toy.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/toy/plot_toy.R">plotting</a>)
* Bayesian Lasso model (<a href="https://github.com/berild/inla-mc/blob/master/lasso/lasso.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/lasso/plot_lasso.R">plotting</a>)
* Spatial Autocorrelation Combined (SAC) model (<a href="https://github.com/berild/inla-mc/blob/master/sem/sem.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/sem/plot_sem.R">plotting</a>)
* Zero-inflated Poisson model (<a href="https://github.com/berild/inla-mc/blob/master/zip/zip.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/zip/output.R">diagnostics</a>)
* Poisson mixture model (<a href="https://github.com/berild/inla-mc/blob/master/pois_mix/pois_mix.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/pois_mix/output.R">diagnostics</a>)

#### Supplementary
* Imputation of missing covariates (<a href="https://github.com/berild/inla-mc/blob/master/missing/missing.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/missing/plot_missing.R">plotting</a>)
* Model-aware Parametric Quantile Regression (<a href="https://github.com/berild/inla-mc/blob/master/pqr/pqr.R">script</a>, <a href="https://github.com/berild/inla-mc/blob/master/pqr/plot_pqr.R">diagnostics</a>)

Furthermore, as some of the examples are time consuming to run, the results from all our examples are available in the <a href="https://github.com/berild/inla-mc/blob/master/sims/">sims</a> folder and the respective figures are available in the <a href="https://github.com/berild/inla-mc/blob/master/figures/">figures</a> folder.

## Genral Usage
To fit a conditional latent Gaussian model with our implementation of the INLA within Monte Carlo methods, a specific function call is required in the respective methods. This call includes the data (`data`), and initial parameters of the proposal distribution (`init = list(x = ... , sigma = ...)`) which could be a mean and covariance matrix but can be anything. The methods also requires some functions as input. The first is the prior (`prior()`) of the parameters being conditioned on $z_c$
```r
prior <- function(x, log = TRUE){
  return(sum(dnorm(x, mean = 0, sd= sqrt(1/.001), log = log)))
}
```
where `log` needs to be a input such that the function can return the log-probability if specified (for in MCMC). Next it requires the proposal distribution
```r
rprop <- function(theta){
  return(mvtnorm::rmvnorm(1, mean=theta[[1]], sigma = theta[[2]]))
}
```
Note, that theta is a self-specified list of parameters controlling the proposal distribution and it is constructed by the `init` parameter in the `inlaIS()` or `inlaAMIS()`. Furthermore, the default adaptation method is mode matching so the algorithm calculates the modes in the <a href="https://github.com/berild/inla-mc/blob/master/inlaMC/genFuncs.R">`calc.theta()`</a> function and the list `theta` will thus contain these modes. If the user wants something else, then simply define your own function called `calc.theta()` similar to what has been done in the <a href="https://github.com/berild/inla-mc/blob/master/pois_mix/pois_mix.R">poisson-mixture example</a>. 
For evaluation, simply add a evaluation point `y` in the function call with the same proposal distribution
```r
dprop <- function(y, theta, log = TRUE){
  return(mvtnorm::dmvnorm(y, mean = theta[[1]], sigma = theta[[2]], log = log))
}
```
Here, we have used the multivariate Gaussian proposal distribution but any proposal distribution can be used as previously mentioned. In MCMC-INLA, the `theta` parameters only specifies the initial state and the hyperparameter of the proposal distribution.
Lastly, a function that uses INLA to fit the conditional model given $z_c$ is required. 

```r
fit.inla <- function(data, beta){
  data$oset = data$x %*% t(beta)
  res = inla(y ~ 1 + offset(oset), data = data)
  return(list(mlik = res$mlik,
               dists = list(intercept = res$marginals.fixed[[1]],
                            precision = res$marginals.hyperpar[[1]])))
}
```
Here, the `beta` parameter represents $z_c$. This function NEEDS to return the conditional posterior marginals of $z_{-c}$ in a list (`dists`) and the conditional marginal likelihood (`mlik`). The inputs are common for all methods (IS-INLA, AMIS-INLA, MCMC-INLA) and needs to contain the `data` first and then the $z_c$ parameter.

To use the IS-INLA method to fit your conditional LGM simply run the following function:
```r
res = inlaIS(data, init, prior, dprop, rprop, fit.inla, N_0, N, ncores)
```
with the inputs described above. Additionally, the input `N_0` specifies the number of samples used in a initial search for a better proposal distribution (not required), `N` is the number of samples of $z_c$ used in the resulting approximation , and `ncores` is the number of CPU cores used.

The AMIS-INLA method is very similar to the IS-INLA function call:
```r
res = inlaAMIS(data, init, prior, dprop, rprop, fit.inla, N_t, N_0, ncores)
```
Here, `N_t` is a vector e.q. `seq(250,500,10)` specifying the number of samples generated by each proposal distribution or prior to each adaptation, and `N_0` is the initial number of samples. 

Lastly, to run our implementation of the MCMC-INLA algorithm simply run the following
```r
res = inlaMH(data, init, prior, dprop, rprop, fit.inla, n.samples, n.burnin, n.thin)
```
where `n.samples` is the number of generated samples, `n.burnin` is the number of samples in the burn-in, and `n.thin` is the number samples in the thinning of the Markov Chain. A implementation of this algorithm can also be found in the `INLABMA::INLAMH()`.

If there is any questions about implementation please send me an email martin.o.berild@ntnu.no.

