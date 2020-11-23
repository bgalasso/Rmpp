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

kmpp <- function(xw, ...)
    UseMethod("kmpp")

kmpp.default <- function(xw, ...) {
    ## Input validation
    stopifnot(is.list(xw))
    ## Disassemble inputs and setup grid
    n <- length(xw)
    aux <- density(xw[[1]], cut = 5, ...)
    grid <- aux$x
    lgrid <- length(grid)
    m <- min(aux$x)
    M <- max(aux$x)
    ## Initialize variables
    y <- Lambda <- matrix(NA, nrow = lgrid, ncol = n)
    bw <- numeric(n)
    ## Step 1 (Learn random measures of warped data)
    cat('Learning...')
    for(i in 1:n) {
        ## cat('\nPoint process ', i)
        aux <- density(xw[[i]], from = m, to = M, ...)
        y[, i] <- aux$y
        bw[i] <- aux$bw
        Lambda[, i] <- cumsum(y[, i] * c(0, diff(grid)))
    }
    ## Steps 2, 3, and 4 (Learn Frechet mean, warp maps, and registration)
    work <- workhorse(xw, grid, Lambda, n, lgrid)
    xreg <- work$xreg
    warp <- work$warp
    F <- work$F
    cat('\nDONE \n')
    ## Organize and return outputs
    class(xreg) <- "mpp"
    class(xw) <- "mpp"
    dens <- list(grid = grid, y = y)
    class(dens) <- c("rcurve", "dens")
    cdf <- list(grid = grid, Lambda = Lambda)
    class(cdf) <- c("rcurve", "cdf")
    frechet <- list(grid = grid, F = F, Lambda = Lambda)
    class(frechet) <- c("rcurve", "fmean")
    warp <- list(grid = grid, warp = warp)
    class(warp) <- c("rcurve", "warp")
    output <- list(xreg = xreg, xw = xw, warp = warp, frechet = frechet,
                   dens = dens, cdf = cdf, bw = bw)
    class(output) <- "kmpp"
  return(output)
}

plot.kmpp <- function(x, type = 'dens' , choose = "all", main = "", bands = TRUE, ...) {
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
    data_pp <- data.frame()
    for(k in 1:n){
      aux <- data.frame(Mean = xw[[k]],
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
    }  else if(is.numeric(choose)) {
      stopifnot(choose > 0 & choose <= n & choose == floor(choose))
      index <- which(data_pp$PProcess == choose)
      data_pp_choose <- data_pp[index,]

      chart <- ggplot() +
        geom_point(data = data_pp_choose,
                   aes(x = Mean, y = PProcess,
                       color =colorRampPalette(c("blue",
                                                 "violet", "red",
                                                 "orange"))(n)[choose]))

      chart <- chart +
        labs(x = "x", y = "Point process") +
        scale_y_continuous(breaks = 1:n) +
        ggtitle(main) +
        theme(legend.position = "none")

      plot(chart)
    }
  } else if(type == 'fmean'){
    cat('Plotting...', '\n')
    x <- x$frechet
    chart <- ggplot(data = data.frame(x = x$grid, y = x$F)) +
      geom_line(aes(x, y), alpha = 1) +
      labs(x = "t", y = "Fr\u{e9}chet mean") +
      ggtitle(main)
    n <- dim(x$Lambda)[2]
    for(i in 1:n) {
      chart <- chart +
        geom_line(data = data.frame(x = x$grid, y = x$Lambda[, i]),
                  aes(x, y), color = i, alpha = 0.1)
    }
    plot(chart)
    cat('DONE', '\n')
  } else if(type == 'dens'){
    cat('Plotting...', '\n')
    x <- x$dens
    n <- dim(x$y)[2]
    if(choose == "all") {
      chart <- ggplot()
      for(i in 1:n) {
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = x$y[, i]),
                    aes(x, y), color = i, alpha = 1)
      }
      chart <- chart +
        labs(x = "x", y = "Density") +
        ggtitle(main)
      plot(chart)
    }
    if(is.numeric(choose)) {
      stopifnot(choose > 0 & choose <= n & choose == floor(choose))
      xw <- x$xw[[choose]]
      chart <- ggplot() +
        geom_area(data = data.frame(x = x$grid, y = x$y[, choose]),
                  aes(x, y), color = choose, fill = choose,
                  alpha = .3) +
        labs(x = "x", y = "Density") +
        ggtitle(main) +
        geom_rug(data = as.data.frame(xw), aes(x = xw), sides = 'b')
      plot(chart)
    }
    cat('DONE', '\n')
  } else if(type == 'cdf') {
    cat('Plotting...', '\n')
    x <- x$cdf
    n <- dim(x$Lambda)[2]
    if(choose == "all") {
      chart <- ggplot()
      for(i in 1:n) {
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = x$Lambda[, i]),
                    aes(x, y), color = i, alpha = 1)
      }
      chart <- chart +
        labs(x = "x", y = "Distribution function") +
        ggtitle(main)
      plot(chart)
    }
    if(is.numeric(choose)) {
      stopifnot(choose > 0 & choose <= n & choose == floor(choose))
      chart <- ggplot() +
        geom_line(data = data.frame(x = x$grid, y = x$Lambda[, choose]),
                  aes(x, y), color = "darkblue", alpha = 1) +
        labs(x = "x", y = "Distribution function") +
        ggtitle(main)
      plot(chart)
    }
    cat('DONE', '\n')
  } else if(type == "warp") {
    cat('Plotting...', '\n')
    x <- x$warp
    n <- dim(x$warp)[2]
    if(choose == "all") {
      chart <- ggplot()
      for(i in 1:n) {
        chart <- chart +
          geom_line(data = data.frame(x = x$grid, y = x$warp[, i]),
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
      xw <- x$xw[[choose]]
      xreg <- x$xreg[[choose]]
      chart <- ggplot() +
        geom_line(data = data.frame(x = x$grid, y = x$warp[, choose]),
                  aes(x, y), color = "darkblue", alpha = 1) +
        labs(x = "x", y = "Warp function") +
        ggtitle(main) +
        geom_abline(slope = 1, intercept = 0,
                    alpha = .1, size = 2) +
        geom_point(data = data.frame(line = rep(0,length(xw)), xw = xw) ,
                   aes(x = xw, y = line), color = 'gray') +
        geom_point(data = data.frame(line = rep(max(x$grid),length(xw)),
                                     xreg = xreg), aes(x = xreg, y = line),
                   color = 'gray')
      plot(chart)
    }
    cat('DONE', '\n')
  }
}
