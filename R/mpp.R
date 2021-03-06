mpp <- function(xw)
    UseMethod("mpp")

plot.mpp <- function(x, choose = "all", main = "", bands = FALSE, ...) {
    ## Input validation
    if(is.character(choose) & choose != "all")
        stop("choose can only be equal to 'all' or to a positive integer.")
    ## Plot
    cat('Plotting...', '\n')
    if(inherits(x,'mpp')){
      n <- length(x)
      if(choose == "all") {
        chart <- ggplot()
        for(i in 1:n) {
          chart <- chart +
            geom_point(data = data.frame(x = x[[i]], y = i),
                       aes(x, y), color = i, alpha = .3)
        }
        chart <- chart +
          labs(x = "x", y = "Point process") +
          scale_y_continuous(breaks = 1:n) +
          ggtitle(main)
        plot(chart)
      }
      else if(is.numeric(choose)) {
        stopifnot(choose > 0 & choose <= n & choose == floor(choose))
        chart <- ggplot() +
          geom_point(data = data.frame(x = x[[choose]], y = choose),
                     aes(x, y), color = choose, alpha = .3) +
          labs(x = "x", y = "Point process") +
          scale_y_continuous(breaks = choose) +
          ggtitle(main)
        plot(chart)
      }
    }

    if(inherits(x, 'lmpp')) {
      n <- length(x)
      if(bands){
        data_pp <- data.frame()
        for(k in 1:n){
          aux <- data.frame(Mean = apply(x[[k]],1,mean),
                            xmin = apply(x[[k]],1,quantile,probs = 0.025,
                                         na.rm = TRUE),
                            xmax = apply(x[[k]],1,quantile,probs = 0.975,
                                         na.rm = TRUE),
                            PProcess = k
          )
          data_pp <- rbind(data_pp,aux)
        }
        if(choose == "all") {
          chart <- ggplot() +
            geom_errorbarh(data = data_pp, aes(y = PProcess,
                                               xmin = xmin, xmax = xmax,
                                               color = PProcess),
                           alpha = .2) +
              geom_point(data = data_pp,
                         aes(x = Mean, y = PProcess, color = PProcess))

          chart <- chart +
            labs(x = "x", y = "Point process") +
            scale_y_continuous(breaks = 1:n) +
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
                           alpha = .1) +
            geom_point(data = data_pp_choose,
                       aes(x = Mean, y = PProcess, color = PProcess))

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
          aux <- data.frame(Mean = apply(x[[k]],1,mean),
                            PProcess = k
          )
          data_pp <- rbind(data_pp,aux)
        }
        if(choose == "all") {
          chart <- ggplot() +
            geom_point(data = data_pp,
                       aes(x = Mean, y = PProcess, color = PProcess))

          chart <- chart +
            labs(x = "x", y = "Point process") +
            scale_y_continuous(breaks = 1:n) +
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
                       aes(x = Mean, y = PProcess, color = PProcess))

          chart <- chart +
            labs(x = "x", y = "Point process") +
            scale_y_continuous(breaks = 1:n) +
            ggtitle(main) +
            theme(legend.position = "none")

          plot(chart)
        }
      }
    }
    cat('DONE', '\n')
}
