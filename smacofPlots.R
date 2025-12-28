
smacofShepardPlot <-
  function(h,
           main = "ShepardPlot",
           fitlines = TRUE,
           colline = "RED",
           colpoint = "BLUE",
           resolution = 100,
           lwd = 2,
           cex = 1,
           pch = 16) {
    maxDelta <- max(h$delta)
    minDelta <- min(h$delta)
    x <- h$delta
    y <- h$dhat
    z <- h$confdist
    plot(
      rbind(cbind(x, z), cbind(x, y)),
      xlim = c(minDelta, maxDelta),
      ylim = c(0, max(c(h$confdist, h$dhat))),
      xlab = "delta",
      ylab = "dhat and dist",
      main = main,
      type = "n"
    )
    points(x,
           z,
           col = colpoint,
           cex = cex,
           pch = pch)
    points(x,
           y,
           col = colline,
           cex = cex,
           pch = pch)
    if (fitlines) {
      for (i in 1:length(x)) {
        lines(x = c(x[i], x[i]), y = c(y[i], z[i]))
      }
    }
    lines(x,
          y,
          type = "l",
          lwd = lwd,
          col = colline)
  }

smacofConfigurationPlot <-
  function(h,
           main = "ConfigurationPlot",
           labels = NULL,
           dim1 = 1,
           dim2 = 2,
           pch = 16,
           col = "RED",
           cex = 1.0) {
    xnew <- h$conf
    if (is.null(labels)) {
      plot(
        xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        pch = pch,
        col = col,
        cex = cex
      )
    }
    else {
      plot(
        xnew[, c(dim1, dim2)],
        xlab = paste("dimension", dim1),
        ylab = paste("dimension", dim2),
        main = main,
        type = "n"
      )
      text(xnew[, c(dim1, dim2)], labels, col = col, cex = cex)
    }
  }

smacofDistDhatPlot <- function(h,
                               fitlines = TRUE,
                               colline = "RED",
                               colpoint = "BLUE",
                               main = "Dist-Dhat Plot",
                               cex = 1.0,
                               lwd = 2,
                               pch = 16) {
  uppe <- max(c(h$confdist, h$dhat))
  plot(
    h$confdist,
    h$dhat,
    xlab = "distance",
    ylab = "disparity",
    xlim = c(0, uppe),
    ylim = c(0, uppe),
    main = main,
    cex = cex,
    pch = pch,
    col = colpoint
  )
  abline(0, 1, col = colline, lwd = lwd)
  if (fitlines) {
    m <- length(h$confdist)
    for (i in 1:m) {
      x <- h$confdist[i]
      y <- h$dhat[i]
      z <- (x + y) / 2
      a <- matrix(c(x, z, y, z), 2, 2)
      lines(a, lwd = lwd)
    }
  }
}

smacofResidualPlot <- function(h, 
                               main = "ResidualPlot",
                               probs = seq(0, 1, 0.25),
                               dats = 1,
                               type = "ecdf",
                               qlines = TRUE,
                               colpoints = "RED",
                               collines = "BLUE",
                               lwdpoints = 2,
                               lwdlines = 2,
                               cex = 1,
                               pch = 16) {
  res <- h$confdist - h$dhat
  res <- switch(dats, res, abs(res), res^2)
  q <- quantile(res, probs)
  n <- length(q)
  e <- switch(type, 
         "ecdf" = ecdf(res),
         "density" = density(res)
         )
  plot(e, col = colpoints, main = main, lwd = lwdpoints)
  for (i in 1:n) {
    abline(v = q[i], col = collines, lwd = lwdlines)
  }
}