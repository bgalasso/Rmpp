# Rmmp guide 

The Rmpp package implements Bayesian semiparametric inferences for modelling, fitting, and registering multiple point processes subject to phase-variation. The main function in the package is the `bmpp` function which is given by  
```
bmpp(xw, support = 1, grid = seq(0, 1, length = 2^8), prior, mcmc, parallel = TRUE)
```
This function implements the model developed in Galasso et al (2020+, arXiv:1812.09607). Some comments on how to use this function are in order. First, multiple point process data should be given in a list whose entries correspond to a single point process.   Here is an example with of data simulated according to this structure: 
```
# Simulate the original point process
n <- 3
l <- 150
m <- rpois(n, l)
x <- list()
for(i in 1:n)
    x[[i]] <- rnorm(m[i], mean = 1 / 2, sd = .15)
# Simulate the warp maps
warp <- rwarp(n)
# Get the warped point process
xw <- list()
for(i in 1:n)
    xw[[i]] <- warp$T(x[[i]],i)
```
Thus, `xw` is our multiple point process in a list form. Next, we need to set the prior and mcmc parameters; following the same specifications as in Galasso et al (2020+), here we consider:  
```
mcmc <- list(nburn = 500, nsave = 4500, nskip = 0, ndisplay = 100)
prior <- list(aa0 = 2, ab0 = 2, kmax = 1000, a0 = 1, b0 = 1)
```
The function implements parallel computing and by default (with `parallel = TRUE`) uses `#cores - 1` cores.   If  `parallel = FALSE`, then the algorithm runs using a single core.
To run our method simply use
```
fit <- bmpp(xw = xw, prior = prior, mcmc = mcmc)
```
The output is stored in the object `fit` that contains the estimated warp maps, the registered point processes, and more.

To visualise the output, one may use the `plot` function which has a specific parameter `type` when the object to plot is of the class `bmpp`; the following options can be used to obtain different targets (respectively Fréchet mean, distribution function, warp maps, and multiple point process):
```
plot(fit, type = 'fmean')
plot(fit, type = 'cdf')       ## try also: plot(fit, type = 'cdf', choose = 1)
plot(fit, type = 'warp')      ## try also: plot(fit, type = 'warp', choose = 1)
plot(fit, type = 'mpp')       ## try also: plot(fit, type = 'mpp', choose = 1)
```
