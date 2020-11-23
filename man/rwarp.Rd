\encoding{UTF-8}
\name{rwarp}
\alias{rwarp}
\alias{rwarp.default}
\title{Random Generation of Warp Maps}
\description{
  Random generation of warp maps.
}
\usage{
## Default S3 method:
rwarp(n, type = "wang-gasser", 
      params = list(lambda = 1, J = 2, alpha = 1, epsilon = 0.05))
}
\arguments{
  \item{n}{number of warp maps to be simulated.}  \item{type}{a
    character string giving the method to be used. This must match one of
    the following: \code{"bend"}, \code{"dp"} or \code{"wang-gasser"}; by default
    \code{type = "wang-gasser"}.}
\item{params}{list with parameters for \code{dp} and \code{wang-gasser}
  type; see details.}
}
\value{
  \item{T}{function with arguments \eqn{(t, i)} giving the family of
    simulated warp functions; the family is indexed \eqn{i}, for
    \eqn{i = 1, \dots, n}.}
  \item{n}{number of simulated warp functions.}
}
\details{
  This function simulates warp functions. If \code{type = "bend"}, then
  trajectories of warp maps are generated as in Galasso \emph{et al.}
  (2020, Section 3.1). For \code{type = "dp"}, the trajectories of the
  warp maps are generated according to a Dirichlet process, centred on a
  uniform distribution, with precision parameter \code{alpha} and
  with approximation parameter \code{epsilon} (by default, \code{params =
  list(alpha = 1, epsilon = 0.05)}). \eqn{\varepsilon}{\epsilon}-DP warp
  maps are constructed by resorting to Muliere and Petrone (1993), that is 
  \deqn{T = \sum_{h=1}^{N_{\varepsilon}} \pi_h \delta_{X_h}, \quad 
    X_h \sim \mbox{Unif[0, 1]}, \quad
    N_{\varepsilon} \sim {\mbox{Poi}(-\alpha \log \varepsilon) + 1},}{
    F = \sum(1,N_\epsilon) (\pi[h] \delta_X[h]),     X[h] ~
    Unif[0, 1],     N_\epsilon ~ Poi(-\alpha log \epsilon) + 1,}
  where the \eqn{\pi_h} are stick-breaking weights. For \code{type =
  "wang-gasser"} trajectories are simulated according to Panaretos and
  Zemel (2016, Section 8.1) (default parameters \code{params =
  list(lambda = 1, J = 2)}).

  The \code{plot} method
  depicts the trajectories of the simulated random warp maps. If
  \code{identity = FALSE}, the reference identity warp map is
  removed from the chart.
}
\references{
  Galasso, B., Zemel, Y. and de Carvalho, M. (2020) Bayesian
  semiparametric modelling of phase-varying point processes. arXiv:
  1812.09607.

  Muliere, P. and Petrone, S. (1993) A Bayesian predictive approach to
  sequential search for an optimal dose: Parametric and nonparametric
  models. \emph{Journal of the Italian Statistical Society}, \bold{2},
  349--364.
  
  Panaretos, V. M. and Zemel, Y. (2016) Amplitude and phase variation of
  point processes. \emph{The Annals of Statistics}, \bold{44}, 771--812.
}
\examples{
## EXAMPLE 1: SIMULATING WARP MAPS
# 1) Bending identity
warp <- rwarp(3, type = "bend")
# 2) Wang--Gasser
plot(warp, main = "Bending identity")
warp <- rwarp(100, type = "wang-gasser")
plot(warp, main = "Wang-Gasser")
# 3) Dirichlet-Uniform  
warp <- rwarp(100, type = "dp", params = list(alpha = 100, epsilon = 0.05))
plot(warp, main = "Dirichlet process")

## EXAMPLE 2: SIMULATING WARPED DATA
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
}
\author{Bastian Galasso-Diaz, Yoav Zemel, and Miguel de Carvalho.}
