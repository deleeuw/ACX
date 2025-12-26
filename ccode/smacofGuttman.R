library(RSpectra)

smacofGuttman <- function(theData, ndim = 2) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  delta <- theData$delta
  wght <- theData$weights
  wsum <- sum(wght)
  dhat <- (theData$delta)^2
  dhat <- dhat * sqrt(wsum / sum(wght * dhat^2))
  bmat <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    bmat[i, j] <- bmat[j, i] <- wght[k] * dhat[k]
  }
  bmat <- -bmat
  diag(bmat) <- -rowSums(bmat)
  ev <- eigs_sym(bmat, ndim, which = "LA")
  evev <- diag(sqrt(pmax(0, ev$values)))
  x <- ev$vectors[, 1:ndim] %*% evev
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    edis[k] <- sum((x[i,] - x[j, ])^2)
  }
  sdd <- sum(wght * edis^2)
  sde <- sum(wght * dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  x <- lbd * x
  sstress <- sum(wght * (dhat - edis)^2) / wsum
  result <- list(
    delta = delta ^ 2,
    dhat = dhat,
    confdist = edis,
    conf = x,
    weightmat = theData$weights,
    sstress = sstress,
    ndim = ndim,
    nobj = nobj,
    iind = theData$iind,
    jind = theData$jind
  )
  class(result) <- c("smacofSSResult", "smacofGuttmanResult")
  return(result)
}
