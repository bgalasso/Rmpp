\encoding{UTF-8}
\name{kmpp}
\alias{kmpp}
\alias{kmpp.default}
\title{Kernel Modeling of Phase-Varying Point Processes}
\description{
  Frequentist nonparametric inference via kernel smoothing
  for modeling, fitting, and registering multiple point processes
  subject to phase-variation.
}
\usage{
## Default S3 method:
kmpp(xw, ...)
}

\arguments{
  \item{xw}{a list with data from multiple point processes from which the
    alignment is to be computed.}
  %% \item{bw}{bandwidth.}
  \item{...}{further arguments for \code{density} methods, such as
    \code{bw} (bandwidth).}
}
\value{
  \item{xreg}{\code{mpp} object containing the registered multiple point
    processes.}
  \item{xw}{\code{mpp} object containing raw warped data.}
  \item{warp}{\code{rcurve} object containing the trajectories of
    the fitted warp maps.}
  \item{frechet}{\code{rcurve} object containing the fitted
    the Frechet mean.}
  \item{dens}{\code{rcurve} object with the \code{n} coordinates of the
    points where the density of the warped data are estimated, along
    with the corresponding estimating density values.}
  \item{cdf}{\code{rcurve} object with the \code{n} coordinates of the
    points where the density of the warped data are estimated, along
    with the corresponding estimating distribution function values.}
  \item{bw}{vector of bandwiths used.}
}
\details{
  This function fits a kernel-based method for for modeling registration
  of multiple point processes as proposed in as proposed in Panaretos
  and Zemel (2016). The method entails modeling the mean measures
  of the phase-varying point processes with a kernel density estimator.
  The \code{plot} method depicts the outputs.
  That is, if \code{type = 'dens'}, \code{'cdf'}, \code{'fmean'},
  \code{'warp'}, or \code{'mpp'} then the densities, distribution
  functions, Fréchet means, warp maps, or the registered point processes
  is respectively depicted. If \code{choose = "all"}, then all outputs
  are depicted; to depict a single object (e.g. point process 1)
  \code{choose} should be set accordingly (e.g. \code{choose = 1}).
}
\references{
  Panaretos, V. M. and Zemel, Y. (2016) Amplitude and phase variation of
  point processes. \emph{The Annals of Statistics}, \bold{44}, 771-812.
}
\examples{
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
fit <- kmpp(xw)
plot(fit, type = 'fmean')
plot(fit, type = 'cdf')       ## try also: plot(fit, type = 'cdf', choose = 1)
plot(fit, type = 'warp')      ## try also: plot(fit, type = 'warp', choose = 1)
plot(fit, type = 'mpp')       ## try also: plot(fit, type = 'mpp', choose = 1)


\dontrun{
## EXAMPLE WITH SIMULATED DATA
# 1) Simulate raw data
N <- 30
points <- 50
mu <- c(1/4, 3/4)
sigma <- c(.02,.03)
w <- c(.2,.8)
lgrid <- 2^8
xx <- seq(0,1,length.out = 2^8)

x <- list()
Npoints <- rpois(N, points)
for(i in 1:N){
  aux <- sample(c(1,2),Npoints[i],replace = T,prob = w)
  x[[i]] <- rnorm(Npoints[i],mean = mu[aux],sd = sigma[aux])
}
# 2) Simulate warp maps
warp <- rwarp(N, type = 'wang-gasser')
# 3) Simulate warped data
xw <- list()
for(i in 1:N)
  xw[[i]] <- warp$T(x[[i]], i)
# 4) Fit model
fit <- kmpp(xw)

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
## Changing years from April to March, instead of January to December)
Temp_Monthly <- split(Temp, f = 'months')
Temp_Yearly <- list()
for(j in 1:N) {
  aux <- c()
  for(i in 1:12) {
    aux <- c(aux,coredata(Temp_Monthly[[4+(i-1)+(j-1)*12]]))
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
for(i in 1:N){
  N_LT[[i]] <- LThresholds[[i]] / 200
  N_HT[[i]] <- (HThresholds[[i]] - 200) / 165
}
#Fit models
fit_low <- kmpp(N_LT)
fit_high <- kmpp(N_HT)

plot(fit_low, type = 'warp')
plot(fit_high, type = 'warp')
}

}
\author{Bastian Galasso-Diaz, Yoav Zemel, and Miguel de Carvalho.}
