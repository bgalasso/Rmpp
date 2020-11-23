##  ========================================================================  ##
##  B. Galasso-Diaz, Y. Zemel, and M. de Carvalho                      	      ##
##  Copyright (C) 2019                                                        ##
##  ------------------------------------------------------------------------  ##
##  This program is free software; you can redistribute it and/or modify      ##
##  it under the terms of the GNU General Public License as published by      ##
##  the Free Software Foundation; either version 2 of the License, or         ##
##  (at your option) any later version.                                       ##
##                                                                            ##
##  This program is distributed in the hope that it will be useful,           ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of            ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             ##
##  GNU General Public License for more details.                              ##
##                                                                            ##
##  You should have received a copy of the GNU General Public License         ##
##  along with this program; if not, a copy is available at                   ##
##  http://www.r-project.org/Licenses/                                        ##
##  ========================================================================  ##

##  ========================================================================  ##
##    Define objects                                                          ##
##  ------------------------------------------------------------------------  ##

workhorse <- function(xw, grid, cdf, n, lgrid)
    UseMethod("workhorse")

workhorse.default <- function(xw, grid, cdf, n, lgrid) {
    ## Numerical inversion is adapted from quantiledensity.R (spatstat package)
    ## Step 1 - Learn empirical Frechet mean
    probs <- seq(0, 1, by = 0.0001)
    lprobs <- length(probs)
    Fi_inv <- matrix(NA, nrow = lprobs, ncol = n)
    for(i in 1:n)
        cdf[, i] <- cdf[, i] / cdf[lgrid, i]
    for(i in 1:lprobs)
        for(j in 1:n) {
            ii <- min(which(cdf[, j] >= probs[i]))
            if(!is.na(ii) && ii >= 1 && ii <= lgrid)
                Fi_inv[i, j] <- grid[ii]
        }
    for(j in 1:n)
      Fi_inv[lprobs,j] = 1
    Finv <- rowMeans(Fi_inv)
    Finv <- Finv / Finv[lprobs] * grid[lgrid]
    F <- c()
    for(i in 1:lgrid) {
        ii <- min(which(Finv >= grid[i]))
        if(!is.na(ii) && ii >= 1 && ii <= lprobs)
            F[i] <- probs[ii]
    }
    F[lgrid] = 1
    ## Step 2 - Learn warp maps
    warp <- matrix(NA, nrow = lgrid, ncol = n)
    for(i in 1:lgrid)
        for(j in 1:n) {
            ii <- min(which(cdf[, j] >= F[i]))
            if(!is.na(ii) && ii >= 1 && ii <= lgrid)
                warp[i, j] <- grid[ii]
        }
    for(j in 1:n)
      warp[lgrid,j] = 1
    ## Step 3 - register points
    xreg <- list()
    for(i in 1:n){
        xr <- c()
        for(j in 1:length(xw[[i]])){
            ii <- min(which(warp[,i] >= xw[[i]][j]))
            if(!is.na(ii) && ii >= 1 && ii <= lgrid)
            xr[j] <- grid[ii]
        }
        xreg[[i]] <- xr
    }
    output <- list(xreg = xreg, warp = warp,Fi_inv = Fi_inv, F = F)
    return(output)
}
