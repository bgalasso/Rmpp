\encoding{UTF-8}
\name{bmpp}
\alias{bmpp}
\alias{bmpp.default}
\title{Bayesian Modelling of Phase-Varying Point Processes}
\description{
  Bayesian semiparametric inference via a Dirichlet--Bernstein prior for
  modelling, fitting, and registering multiple point processes subject
  to phase-variation.
}
\usage{
## Default S3 method:
bmpp(xw, support = 1, grid = seq(0, 1, length = 2^8), prior, mcmc, parallel = TRUE)
}
\arguments{
  \item{xw}{a list with data from multiple point processes from which the
    alignment is to be computed.}
   \item{support}{an integer number giving the support of the random
     density, \code{1 = [0,1]}, \code{2 = (0, +Inf]}, and
     \code{3 = (-Inf,+Inf)}. Depending on this,
     the data is transformed to lie in the unit interval.}
   %% \item{ngrid}{number of grid points where the density estimate is
   %%   evaluated. This is only used if dimension of \code{y}
   %%   is lower or equal than 2. The default value is 1000.}
   \item{grid}{vector of grid points where the density estimate is
     evaluated. The default value is NULL and the grid is chosen
     according to the range of the data.}
   \item{prior}{a list with prior information for the Dirichlet--Bernstein
     polynomial prior. The list includes the following parameters:
     \code{aa0} and \code{ab0} giving the hyperparameters for
     prior distribution of the precision parameter of the Dirichlet process
     prior, \code{alpha} giving the value of the precision parameter (it
     must be specified if \code{aa0} is missing, see details
     below), \code{a0} and \code{b0} giving the parameters of the
     beta centering distribution of the DP prior, and
     \code{kmax} giving the maximum value of the discrete uniform
     prior for the degree of the Bernstein polynomial.
   }
  \item{mcmc}{a list giving the MCMC parameters. The list must include
    the following integers: \code{nburn} giving the number of burn-in
    scans, \code{nskip} giving the thinning interval, \code{nsave} giving
    the total number of scans to be saved, and \code{ndisplay} giving
    the number of saved scans to be displayed on screen (the function reports
    on the screen when every \code{ndisplay} iterations have been carried
    out).}
  \item{parallel}{logical value indicating whether single core or multiple cores are used.  The default value is \code{TRUE}.}
  %% \item{grid}{grid on which the posterior object of interest is to be
  %%   evaluated; by default \code{grid = seq(0, 1, length = 2^8)}.}
  %% \item{parallel}{logical; if \code{TRUE}, the alignment process will be
  %%   conducted using \eqn{n - 1} cores, where \eqn{n} is the total number
  %%   of (virtual) cores.}
}

\value{
  \item{xreg}{list containing posterior simulated trajectories of
    the registered point process.}
  \item{xreg.hat}{list containing posterior mean of
    the registered point process.}
  \item{xw}{\code{mpp} object containing raw warped data.}
  \item{warp}{3D-array (of size \code{ngrid} \eqn{* n *}
    \code{nsave}) containing posterior simulated trajectories of the warp
    maps, where \eqn{n} is the number of point processes.}
  \item{frechet}{matrix containing posterior simulated trajectories of
    the Frechet mean.}
  \item{x}{list with raw data.}
}
\details{
  This function fits a Bayesian semiparametric method for modeling
  registration of multiple point processes as proposed in Galasso,
  Zemel, and de Carvalho (2020). The method entails modeling the mean measures
  of the phase-varying point processes with a Dirichlet--Bernstein 
  prior, which induces a prior in the space of all warp functions.
  The \code{plot} method depicts the outputs.
  That is, if \code{type = 'dens'}, \code{'cdf'}, \code{'fmean'},
  \code{'warp'}, or \code{'mpp'} then the densities, distribution
  functions, Fréchet means, warp maps, or the registered point processes
  is respectively depicted. If \code{choose = "all"}, then all outputs
  are depicted; to depict a single object (e.g. point process 1)
  \code{choose} should be set accordingly (e.g. \code{choose = 1}).
  If \code{bands = TRUE} the bands are depicted; if \code{last = TRUE}
  only the last 100 trajectories are represented (the latter option is
  not available for \code{mpp}).
}
\references{
  Galasso, B., Zemel, Y. and de Carvalho, M. (2020) Bayesian
  semiparametric modelling of phase-varying point processes. arXiv:
  1812.09607.
}
\examples{
\dontrun{
## EXAMPLE WITH SIMULATED DATA
# 1) Simulate raw data
n <- 3
l <- 150
m <- rpois(n, l)
x <- list()
for(i in 1:n)
    x[[i]] <- rnorm(m[i], mean = 1 / 2, sd = .15)
# 2) Simulate warp maps
warp <- rwarp(n)
# 3) Simulate warped data
xw <- list()
for(i in 1:n)
    xw[[i]] <- warp$T(x[[i]], i)
# 4) Fit model

fit <- bmpp(xw)

plot(fit, type = 'fmean')
plot(fit, type = 'cdf')       ## try also: plot(fit, type = 'cdf', choose = 1)
plot(fit, type = 'warp')      ## try also: plot(fit, type = 'warp', choose = 1)
plot(fit, type = 'mpp')       ## try also: plot(fit, type = 'mpp', choose = 1)
}
\dontrun{
## EXAMPLE WITH SIMULATED DATA
# 1) Simulate raw data
N <- 30
points <- 50
mu <- c(1/4, 3/4)
sigma <- c(.02, .03)
w <- c(.2, .8)
lgrid <- 2^8
xx <- seq(0, 1, length.out = 2^8)

x <- list()
Npoints <- rpois(n, points)
for(i in 1:n){
  aux <- sample(c(1,2), Npoints[i], replace = T, prob = w)
  x[[i]] <- rnorm(Npoints[i], mean = mu[aux], sd = sigma[aux])
}
# 2) Simulate warp maps
warp <- rwarp(n, type = 'wang-gasser')
# 3) Simulate warped data
xw <- list()
for(i in 1:n)
  xw[[i]] <- warp$T(x[[i]], i)
# 4) Fit model

fit <- bmpp(xw)

class(x) <- 'mpp'
plot(x)
plot(fit$xw)
plot(fit, type = 'mpp')
}

\dontrun{
## EXAMPLE WITH REAL DATA
data(santiago)
## Data as a time series
require(xts)
Temp <- xts(santiago[[2]],order.by = santiago[[1]])
N <- floor(length(Temp)/365)
## Changing years from April to March, instead of January to December
Temp_Monthly <- split(Temp, f = 'months')
Temp_Yearly <- list()
for(j in 1:N) {
  aux <- c()
  for(i in 1:12) {
    aux <- c(aux, coredata(Temp_Monthly[[4+(i-1)+(j-1)*12]]))
  }
  Temp_Yearly[[j]] <- aux
}
HThresholds <- LThresholds <-  list()
for(i in 1:N) {
  threshold1 <- quantile(Temp_Yearly[[i]],p = .05)
  threshold2 <- quantile(Temp_Yearly[[i]],p = .95)
  LThresholds[[i]] <- which(Temp_Yearly[[i]] <= threshold1)
  HThresholds[[i]] <- which(Temp_Yearly[[i]] >= threshold2)
}
N_HT <- N_LT <- list()
for(i in 1:N) {
  N_LT[[i]] <- LThresholds[[i]] / 200
  N_HT[[i]] <- (HThresholds[[i]] - 200) / 165
}
#Fit models

fit_low <- bmpp(N_LT)
fit_high <- bmpp(N_HT)

plot(fit_low, type = 'warp')
plot(fit_high, type = 'warp')
}

}
\author{Bastian Galasso-Diaz, Yoav Zemel, and Miguel de Carvalho.}
