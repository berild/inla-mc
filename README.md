# Importance Sampling with the Integrated Nested Laplace Approximation

This repository contains the code used in the Importance Sampling with the Integrated Nested Laplace Approximation paper. The implementation of the Importance Sampling with the Integrated Nested Laplace Approximation (IS-INLA) can be found in <a href="https://github.com/berild/inla-mc/blob/master/inlaIS.R">inlaIS.R</a>, the Adaptive Multiple Importance Sampling with the Integrated Nested Laplace Approximation (AMIS-INLA) in <a href="https://github.com/berild/inla-mc/blob/master/inlaAMIS.R">inlaAMIS.R</a>, and the Markov Chain Monte Carlo with the Integrated Nested Laplace Approximation (MCMC-INLA) in <a href="https://github.com/berild/inla-mc/blob/master/inlaMH.R">inlaMH.R</a>. In addition, general functions used in all algorithms are collected in <a href="https://github.com/berild/inla-mc/blob/master/genFuncs.R">genFuncs.R</a>. 


## Usage

```r
inlaIS(data, init, prior, dprop, rprop, fit.inla, N_0, N, ncores)
```

```r
inlaAMIS(data, init, prior, dprop, rprop, fit.inla, N_t, N_0, ncores)
```

```r
 inlaMH(data, init, prior, dprop, rprop, fit.inla, n.samples, n.burnin, n.thin)
```

```r
prior <- function(x){

}
```

```r
dprop <- function(x){

}
```

```r
rprop <- function(x){

}
```

```r
fit.inla <- function(x){

}
```

### Simple Example
