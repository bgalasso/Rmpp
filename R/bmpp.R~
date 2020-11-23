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

bmpp <- function(xw, support = 1, grid = seq(0,1, length = 2^8), prior, mcmc,
                 parallel = TRUE)
    UseMethod("bmpp")

bmpp.default <- function(xw, support = 1, grid = seq(0, 1, length = 2^8),
                         prior = list(aa0 = 2, ab0 = 2, kmax = 1000, a0 = 1, b0 = 1), 
                         mcmc = list(nburn = 100, nsave = 1000, nskip = 0, ndisplay = 100),
                         parallel = TRUE) {
    ## Input validation and setup grid
    stopifnot(is.list(xw), is.list(prior), is.list(mcmc), is.vector(grid))
    ## Step 1 (Learn random measures of warped data)
    n <- length(xw)
    nsave <- mcmc$nsave
    lgrid <- length(grid)
    CORES <- detectCores() - 1

    y <- array(NA, dim = c(lgrid, nsave, n))
    Lambda <- array(NA, dim = c(lgrid, nsave, n))

    cat('Learning...\n')
    if(parallel){
      bernfun <- function(y){
        aux <- BDPdensity(y, support = support,
                          ngrid = lgrid, prior = prior,
                          mcmc = mcmc,
                          grid = grid,
                          state = NULL,
                          status = TRUE)
        aux <- matrix(aux$fun, nrow = lgrid)
        Lambda.f <- matrix(NA,nrow  = lgrid, ncol = nsave)
        for(t in 1:nsave)
          Lambda.f[, t] <- cumsum(aux[,t] * c(0, diff(grid)))
        res <- list(dens = aux, dist = Lambda.f)
        return(res)
      }
      cat('\nPoint Processes... \n')
      bd.output <- pbmclapply(X = xw, FUN = bernfun,mc.cores = CORES)
      for(i in 1:n){
        y[,,i] <- bd.output[[i]]$dens
        Lambda[,,i] <- bd.output[[i]]$dist
      }
    } else {
      for(i in 1:n) {
          cat('\nPoint process ', i,'\n')
          aux <- BDPdensity(xw[[i]], support = support,
                            ngrid = lgrid, prior = prior,
                            mcmc = mcmc,
                            grid = grid,
                            state = NULL,
                            status = TRUE)
          y[, , i] <- matrix(aux$fun, nrow = lgrid)
          for(t in 1:nsave)
              Lambda[, t, i] <- cumsum(y[, t, i] * c(0, diff(grid)))
      }
    }

    ## Steps 2, 3, and 4 (Learn Frechet mean, warp maps, and registration)
    probs <- seq(0, 1, by = 0.0001)
    lprobs <- length(probs)
    xreg <- vector(mode = "list", length = n)
    warp <- array(NA, dim = c(lgrid,n,nsave))
    F <- matrix(NA, nrow = lgrid, ncol = nsave)
    Fi_inv <- array(NA, dim = c(lprobs,n,nsave))

    if(parallel){
      work.parallel <- function(y){
        res <- workhorse(xw, grid, y, n, lgrid)
        return(res)
      }
      Lambda.list <- list()
      for(t in 1:nsave)
        Lambda.list[[t]] <- Lambda[,t,]

      work.output <- pbmclapply(X =Lambda.list,FUN = work.parallel,
                                mc.cores = CORES)
      for(t in 1:nsave){
        for(k in 1:n)
          xreg[[k]] <- cbind(xreg[[k]],work.output[[t]]$xreg[[k]])
        warp[, , t] <- work.output[[t]]$warp
        Fi_inv[, , t] <- work.output[[t]]$Fi_inv
        F[, t] <- work.output[[t]]$F
      }
    } else {
      pb <- txtProgressBar(min = 0, max = nsave, style = 3)
      for(t in 1:nsave) {
        setTxtProgressBar(pb, t)
        work <- workhorse(xw, grid, Lambda[, t, ], n, lgrid)
        for(k in 1:n)
          xreg[[k]] <- cbind(xreg[[k]],work$xreg[[k]])
        warp[, , t] <- work$warp
        Fi_inv[, , t] <- work$Fi_inv
        F[, t] <- work$F
      }
    }

    xreg.hat <- lapply(xreg,FUN = function(x){apply(x,1,mean)})
    cat('\nDONE \n')
    ## Organize and return outputs
    class(xreg) <- "lmpp"
    class(xreg.hat) <- 'mpp'
    class(xw) <- "mpp"
    dens <- list(grid = grid, y = y)
    class(dens) <- c("rcurve", "bdens")
    cdf <- list(grid = grid, Lambda = Lambda)
    class(cdf) <- c("rcurve", "bcdf")
    frechet <- list(grid = grid, F = F, Lambda = Lambda)
    class(frechet) <- c("rcurve", "bfmean")
    warp <- list(grid = grid, warp = warp)
    class(warp) <- c("rcurve", "bwarp")
    output <- list(xreg = xreg, xreg.hat = xreg.hat ,xw = xw, warp = warp,
                   frechet = frechet, dens = dens, cdf = cdf,
                   Fi_inv = Fi_inv)
    class(output) <- "bmpp"
    return(output)
}

plot.bmpp <- function(x,type = 'dens', choose = "all", main = "",
                      bands = FALSE, last = FALSE, ...){
  ## Input validation
  if(is.character(choose) & choose != "all")
    stop("choose can only be equal to 'all' or to a positive integer.")
  if(!(type == 'fmean' | type == 'dens' | type == 'cdf' | type == 'warp' |
       type == 'mpp'))
    stop("type can only be one of the following: 'mpp','fmean', 'dens', 'cdf' or 'warp'")

  ## Plot
  if(type == "mpp"){
    xw <- x$xreg
    n <- length(xw)
    if(bands){
      data_pp <- data.frame()
      for(k in 1:n){
        aux <- data.frame(Mean = apply(xw[[k]],1,mean),
                          xmin = apply(xw[[k]],1,quantile,probs = 0.025,
                                       na.rm = TRUE),
                          xmax = apply(xw[[k]],1,quantile,probs = 0.975,
                                       na.rm = TRUE),
                          PProcess = k
        )
        data_pp <- rbind(data_pp,aux)
      }
      if(choose == "all") {
        chart <- ggplot() +
          geom_errorbarh(data = data_pp, aes(y = PProcess,
                                             xmin = xmin, xmax = xmax),
                         color = 'gray',
                         alpha = .2) +
          geom_point(data = data_pp,
                     aes(x = Mean, y = PProcess, color = as.factor(PProcess)))

        chart <- chart +
          labs(x = "x", y = "Point process") +
          scale_y_continuous(breaks = 1:n) +
          scale_color_manual(values = colorRampPalette(c("blue", "violet",
                                                         "red", "orange"))(n)) +
          ggtitle(main) +
          theme(legend.position = "none")

        plot(chart)
      }
      else if(is.numeric(choose)) {
        stopifnot(choose > 0 & choose <= n & choose == floor(choose))
        index <- which(data_pp$PProcess == choose)
        data_pp_choose <- data_pp[index,]

        chart <- ggplot() +
          geom_errorbarh(data = data_pp_choose, aes(y = PProcess,
                                                    xmin = xmin, xmax = xmax,
                                                    color = PProcess),
                         alpha = .2) +
          geom_point(data = data_pp_choose,
                     aes(x = Mean, y = PProcess,
                         color = PProcess))

        chart <- chart +
          labs(x = "x", y = "Point process") +
          scale_y_continuous(breaks = 1:n) +
          ggtitle(main) +
          theme(legend.position = "none")

        plot(chart)
      }
    } else {
      data_pp <- data.frame()
      for(k in 1:n){
        aux <- data.frame(Mean = apply(xw[[k]],1,mean),
                          PProcess = k
        )
        data_pp <- rbind(data_pp,aux)
      }
      if(choose == "all") {
        chart <- ggplot() +
          geom_point(data = data_pp,
                     aes(x = Mean, y = PProcess, color = as.factor(PProcess)))

        chart <- chart +
          labs(x = "x", y = "Point process") +
          scale_y_continuous(breaks = 1:n) +
          scale_color_manual(values = colorRampPalette(c("blue", "violet",
                                                         "red", "orange"))(n)) +
          ggtitle(main) +
          theme(legend.position = "none")

        plot(chart)
      }
      else if(is.numeric(choose)) {
        stopifnot(choose > 0 & choose <= n & choose == floor(choose))
        index <- which(data_pp$PProcess == choose)
        data_pp_choose <- data_pp[index,]

        chart <- ggplot() +
          geom_point(data = data_pp_choose,
                     aes(x = Mean, y = PProcess,
                         color = colorRampPalette(c("blue", "violet",
                                                    "red", "orange"))(n)[k]))

        chart <- chart +
          labs(x = "x", y = "Point process") +
          scale_y_continuous(breaks = 1:n) +
          ggtitle(main) +
          theme(legend.position = "none")

        plot(chart)
      }
    }
  } else if(type == "fmean") {
    x <- x$frechet
    F_hat <- apply(x$F,1,mean)
    cat('Plotting...', '\n')
    if(bands){
      Fl <- apply(x$F, 1, quantile, probs = 0.025, na.rm = TRUE)
      Fu <- apply(x$F, 1, quantile, probs = 0.975, na.rm = TRUE)
      credible.bands <- data.frame(X = c(x$grid,rev(x$grid)),
                                   Y = c(Fl,rev(Fu)))

      chart <- ggplot(data = data.frame(x = x$grid, y = F_hat)) +
        geom_line(aes(x, y), alpha = 1) +
        labs(x = "x", y = "Fr\u{e9}chet mean") +
        ggtitle(main) +
        geom_polygon(data = credible.bands,
                     aes(x = X, y = Y, fill = 3),
                     alpha = 0.3) +
        theme(legend.position = "none")
    } else if(last) {
      K <- dim(x$F)[2]
      stopifnot(K >= 100)
      chart <- ggplot()
      for(k in (K-99):K){
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = x$F[,k]),
                    aes(x,y), color = 'gray', alpha = 0.7)
      }
      chart <- chart +
        labs(x = "x", y = "Density") +
        ggtitle(main)
    } else {
      chart <- ggplot(data = data.frame(x = x$grid, y = F_hat)) +
        geom_line(aes(x, y), alpha = 1) +
        labs(x = "x", y = "Fr\u{e9}chet mean") +
        ggtitle(main)
      n <- dim(x$Lambda)[3]
      for(i in 1:n) {
        Lambda <- apply(x$Lambda[,,i],1,mean)
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = Lambda),
                    aes(x, y), color = i, alpha = 0.1)
      }
    }
    plot(chart)
    cat('DONE', '\n')
  } else if(type == "dens") {
    x <- x$dens
    cat('Plotting...', '\n')
    n <- dim(x$y)[3]
    if(choose == "all") {
      chart <- ggplot()
      for(i in 1:n) {
        denb <- apply(x$y[,,i],1,mean)
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = denb),
                    aes(x, y), color = i, alpha = 1)
      }
      chart <- chart +
        labs(x = "x", y = "Density") +
        ggtitle(main)
      plot(chart)
    }
    if(is.numeric(choose)) {
      stopifnot(choose > 0 & choose <= n & choose == floor(choose))
      denb <- apply(x$y[,,choose],1,mean)
      xw <- x$xw[[choose]]
      if(bands){
        densl <- apply(x$y[,,choose], 1, quantile, probs = 0.025, na.rm = TRUE)
        densu <- apply(x$y[,,choose], 1, quantile, probs = 0.975, na.rm = TRUE)
        credible.bands <- data.frame(X = c(x$grid,rev(x$grid)),
                                     Y = c(densl,rev(densu)))
        chart <- ggplot(data = data.frame(x = x$grid, y = denb),
                        aes(x, y)) +
          geom_line() +
          labs(x = "x", y = "Density") +
          ggtitle(main) +
          geom_polygon(data = credible.bands,
                       aes(x = X, y = Y, fill = 3),
                       alpha = 0.3) +
          theme(legend.position = "none") +
          geom_rug(data = as.data.frame(xw), aes(x = xw, y = 0), sides = 'b')
      } else if(last) {
        K <- dim(x$y)[2]
        stopifnot(K >= 100)
        chart <- ggplot()
        for(k in (K-99):K){
          chart <- chart +
            geom_line(data = data.frame(x = x$grid, y = x$y[,k,choose]),
                      aes(x,y), color = 'gray', alpha = 0.7)
        }
        chart <- chart +
          labs(x = "x", y = "Density") +
          ggtitle(main)
      } else {
        chart <- ggplot() +
          geom_area(data = data.frame(x = x$grid, y = denb),
                    aes(x, y), color = choose, fill = choose,
                    alpha = .3) +
          labs(x = "x", y = "Density") +
          ggtitle(main) +
          geom_rug(data = as.data.frame(xw), aes(x = xw, y = 0), sides = 'b')
      }
      plot(chart)
    }
    cat('DONE', '\n')
  } else if(type == "cdf") {
    x <- x$cdf
    cat('Plotting...', '\n')
    n <- dim(x$Lambda)[3]
    if(choose == "all") {
      chart <- ggplot()
      for(i in 1:n) {
        Lambda <- apply(x$Lambda[,,i],1,mean)
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = Lambda),
                    aes(x, y), color = i, alpha = 1)
      }
      chart <- chart +
        labs(x = "x", y = "Distribution function") +
        ggtitle(main)
      plot(chart)
    }
    if(is.numeric(choose)) {
      stopifnot(choose > 0 & choose <= n & choose == floor(choose))
      Lambda <- apply(x$Lambda[,,choose],1,mean)
      if(bands){
        laml <- apply(x$Lambda[,,choose], 1, quantile, probs = 0.025, na.rm = TRUE)
        lamu <- apply(x$Lambda[,,choose], 1, quantile, probs = 0.975, na.rm = TRUE)
        credible.bands <- data.frame(X = c(x$grid,rev(x$grid)),
                                     Y = c(laml,rev(lamu)))
        chart <- ggplot() +
          geom_line(data = data.frame(x = x$grid, y = Lambda),
                    aes(x, y), color = "darkblue", alpha = 1) +
          labs(x = "x", y = "Distribution function") +
          ggtitle(main) +
          geom_polygon(data = credible.bands,
                       aes(x = X, y = Y, fill = 3),
                       alpha = 0.3) +
          theme(legend.position = "none")
      } else if(last) {
        K <- dim(x$Lambda)[2]
        stopifnot(K >= 100)
        chart <- ggplot()
        for(k in (K-99):K){
          chart <- chart +
            geom_line(data = data.frame(x = x$grid, y = x$Lambda[,k,choose]),
                      aes(x,y), color = 'gray', alpha = 0.7)
        }
        chart <- chart +
          labs(x = "x", y = "Density") +
          ggtitle(main)
      } else {
        chart <- ggplot() +
          geom_line(data = data.frame(x = x$grid, y = Lambda),
                    aes(x, y), color = "darkblue", alpha = 1) +
          labs(x = "x", y = "Distribution function") +
          ggtitle(main)
      }
      plot(chart)
    }
    cat('DONE', '\n')
  } else if(type == "warp") {
    x <- x$warp
    cat('Plotting...', '\n')
    n <- dim(x$warp)[2]
    if(choose == "all") {
      chart <- ggplot()
      for(i in 1:n) {
        warp <- apply(x$warp[,i,],1,mean)
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = warp),
                    aes(x, y), color = i, alpha = 1)
      }
      chart <- chart +
        labs(x = "x", y = "Warp function") +
        ggtitle(main) +
        geom_abline(slope = 1, intercept = 0,
                    alpha = .1, size = 2)
      plot(chart)
    }
    if(is.numeric(choose)) {
      stopifnot(choose > 0 & choose <= n & choose == floor(choose))
      warp <- apply(x$warp[,choose,],1,mean)
      xw <- x$xw[[choose]]
      xreg <- x$xreg.hat[[choose]]
      if(bands){
        warpl <- apply(x$warp[,choose,], 1, quantile, probs = 0.025, na.rm = TRUE)
        warpu <- apply(x$warp[,choose,], 1, quantile, probs = 0.975, na.rm = TRUE)
        credible.bands <- data.frame(X = c(x$grid,rev(x$grid)),
                                     Y = c(warpl,rev(warpu)))
        chart <- ggplot() +
          geom_line(data = data.frame(x = x$grid, y = warp),
                    aes(x, y), color = "darkblue", alpha = 1) +
          labs(x = "x", y = "Warp function") +
          ggtitle(main) +
          geom_abline(slope = 1, intercept = 0,
                      alpha = .1, size = 2) +
          geom_polygon(data = credible.bands,
                       aes(x = X, y = Y, fill = 3),
                       alpha = 0.3) +
          theme(legend.position = "none")
        chart <- chart +
          geom_point(data = data.frame(line = rep(0,length(xw)), xw = xw) ,
                     aes(x = xw, y = line), color = 'gray') +
          geom_point(data = data.frame(line = rep(max(x$grid),length(xw)),
                                       xreg = xreg), aes(x = xreg, y = line),
                     color = 'gray')
      } else if(last) {
        K <- dim(x$warp)[3]
        stopifnot(K >= 100)
        chart <- ggplot()
        for(k in (K-99):K){
          chart <- chart +
            geom_line(data = data.frame(x = x$grid, y = x$warp[,choose,k]),
                      aes(x,y), color = 'gray', alpha = 0.7)
        }
        chart <- chart +
          labs(x = "x", y = "Density") +
          ggtitle(main)
      } else {
        chart <- ggplot() +
          geom_line(data = data.frame(x = x$grid, y = warp),
                    aes(x, y), color = "darkblue", alpha = 1) +
          labs(x = "x", y = "Warp function") +
          ggtitle(main) +
          geom_abline(slope = 1, intercept = 0,
                      alpha = .1, size = 2)
        chart <- chart +
          geom_point(data = data.frame(line = rep(0,length(xw)), xw = xw) ,
                     aes(x = xw, y = line), color = 'gray') +
          geom_point(data = data.frame(line = rep(max(x$grid),length(xw)),
                                       xreg = xreg), aes(x = xreg, y = line),
                     color = 'gray')
      }
      plot(chart)
    }
    cat('DONE', '\n')
  }
}



