library(RSpectra)

smacofTorgerson <- function(theData, ndim = 2) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  wght <- theData$weights
  dmat <- matrix(0, nobj, nobj)
  dhat <- (theData$delta)^2
  mdel <- mean(dhat)
  dmat <- mdel * (1 - diag(nobj))
  wmat <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- theData$iind[k]
    j <- theData$jind[k]
    dmat[i, j] <- dmat[j, i] <- dhat[k]
    wmat[i, j] <- wmat[j, i] <- wght[k]
  }
  wsum <- sum(wmat)
  dmat <- dmat * sqrt(wsum / sum(wmat * dmat^2))
  dr <- apply(dmat, 1, mean)
  dm <- mean(dmat)
  bmat <- -(dmat - outer(dr, dr, "+") + dm) / 2
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
    delta = theData$delta ^ 2,
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
  class(result) <- c("smacofSSResult", "smacofTorgersonResult")
  return(result)
  }
