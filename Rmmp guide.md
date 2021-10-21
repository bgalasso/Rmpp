# Rmmp guide 

Rmpp package implement the Bayesian semi-parametric inferences for modelling, fitting and registering multiple point processes subject to phase-variation. The main function in the package is `bmpp` function which is 
```
bmpp(xw, support = 1, grid = seq(0, 1, length = 2^8), prior, mcmc, parallel = TRUE)
```
which implement the model developed in the paper. As a guide, we explain here how use this function. First we need to get our multiple point processes data as a list, which every element in that list correspond to a single point process, for example:
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
now `xw` is our multiple point process in a list form, so then we need to define the prior and mcmc parameters. As in the paper, we use some standard configuration for both
```
mcmc <- list(nburn = 500, nsave = 4500, nskip = 0, ndisplay = 100)
prior <- list(aa0 = 2, ab0 = 2, kmax = 1000, a0 = 1, b0 = 1)
```
and finally, we need to define if we can run this in parallel or not. By default, the algorithm run using `#cores - 1` (`parallel = TRUE`), in the case that `parallel = FALSE`, the the algorithm run in a single-core.
To run our algorithm, you need just execute the following command:
```
fit <- bmpp(xw = xw, prior = prior, mcmc = mcmc)
```
The output was stored in the object `fit` which has as output the estimated warp mpas and registered point processes, among others.

To visualize the output, you could use the `plot` function which has a specific parameter `type` when the object to plot is the class `bmpp`, so you have the folowwing options:
```
plot(fit, type = 'fmean')
plot(fit, type = 'cdf')       ## try also: plot(fit, type = 'cdf', choose = 1)
plot(fit, type = 'warp')      ## try also: plot(fit, type = 'warp', choose = 1)
plot(fit, type = 'mpp')       ## try also: plot(fit, type = 'mpp', choose = 1)
```
