

makeMDSData <- function(delta, weights = NULL) {
  nobj <- attr(delta, "Size")
  if (is.null(weights)) {
    weights <- as.dist(1 - diag(nobj))
  }
  theData <- NULL
  k <- 1
  for (j in 1:(nobj - 1)) {
    for (i in (j + 1):nobj) {
      if ((weights[k] > 0) &&
          (!is.na(weights[k])) && (!is.na(delta[k]))) {
        theData <- rbind(theData, c(i, j, delta[k], 0, weights[k]))
      }
      k <- k + 1
    }
  }
  colnames(theData) <- c("i", "j", "delta", "blocks", "weights")
  ndat <- nrow(theData)
  theData <- theData[order(theData[, 3]), ]
  dvec <- theData[, 3]
  k <- 1
  repeat {
    m <- length(which(dvec == dvec[k]))
    theData[k, 4] <- m
    k <- k + m
    if (k > ndat) {
      break
    }
  }
  result <- list(
    iind = theData[, 1],
    jind = theData[, 2],
    delta = theData[, 3],
    blocks = theData[, 4],
    weights = theData[, 5],
    nobj = nobj,
    ndat = ndat
  )
  class(result) <- "smacofSSData"
  return(result)
}

fromMDSData <- function(theData) {
  ndat <-theData$ndat
  nobj <- theData$nobj
  delta <- matrix(0, nobj, nobj)
  weights <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    delta[i, j] <- delta[j, i] <- theData$delta[k]
    weights[i, j] <- weights[j, i] <- theData$weights[kdelta]
  }
  return(list(delta = as.dist(delta), weights = as.dist(weights)))
}


smacofRandomConfiguration <- function(theData, ndim = 2) {
  nobj <- theData$nobj
  x <- matrix(rnorm(nobj * ndim), nobj, ndim)
  return(apply(x, 2, function(x) x - mean(x)))
}

matrixPrint <- function(x,
                        digits = 6,
                        width = 8,
                        format = "f",
                        flag = "+") {
  print(noquote(
    formatC(
      x,
      digits = digits,
      width = width,
      format = format,
      flag = flag
    )
  ))
}
