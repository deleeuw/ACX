library(RSpectra)

smacofElegant <- function(theData,
                          ndim = 2,
                          itmax = 1000,
                          eps = 1e-10,
                          quick = FALSE,
                          verbose = FALSE) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  wght <- theData$weights
  wsum <- sum(wght)
  dhat <- (theData$delta)^2
  dhat <- dhat * sqrt(wsum / sum(wght * dhat^2))
  cinit <- smacofDoubleCenter(theData)
  lbda <- smacofBound(theData)
  evev <- eigs_sym(cinit, ndim)
  cinit <- tcrossprod(evev$vectors %*% diag(sqrt(evev$values)))
  dold <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    dold[k] = cinit[i, i] + cinit[j, j] - 2 * cinit[i, j]
  }
  sold <- sum(wght * (dhat - dold)^2) / wsum
  itel <- 1
  cold <- cinit
  repeat {
    cmaj <- matrix(0, nobj, nobj)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      cmaj[i, j] <- cmaj[j, i] <- wght[k] * (dhat[k] - dold[k])
    }
    cmaj <- -cmaj
    diag(cmaj) <- -rowSums(cmaj)
    cmaj <- cold + cmaj / lbda
    evev <- eigs_sym(cmaj, ndim)
    xvev <- evev$vectors %*% diag(sqrt(evev$values))
    cnew <- tcrossprod(xvev)
    dnew <- rep(0, ndat)
    for (k in 1:ndat) {
      i <- iind[k]
      j <- jind[k]
      dnew[k] = cnew[i, i] + cnew[j, j] - 2 * cnew[i, j]
    }
    snew <- sum(wght * (dhat - dnew)^2) / wsum
    if (verbose) {
      cat(
        "itel ",
        formatC(itel, format = "d"),
        "sold ",
        formatC(sold, digits = 10, format = "f"),
        "snew ",
        formatC(snew, digits = 10, format = "f"),
        "\n"
      )
    }
    if ((itel == itmax) || ((sold - snew) < eps)) {
      break
    }
    itel <- itel + 1
    dold <- dnew
    sold <- snew
    cold <- cnew
  }
  result <- list(
    delta = theData$delta ^ 2,
    dhat = dhat,
    confdist = dnew,
    conf = xvev,
    weightmat = theData$weights,
    stress = snew,
    ndim = ndim,
    niter = itel,
    nobj = nobj,
    iind = iind,
    jind = jind
  )
  class(result) <- c("smacofSSResult", "smacofElegantResult")
  return(result)
}

smacofDoubleCenter <- function(theData) {
  nobj <- theData$nobj
  ndat <- theData$ndat
  iind <- theData$iind
  jind <- theData$jind
  wght <- theData$weights
  wsum <- sum(wght)
  dhat <- (theData$delta)^2
  dhat <- dhat * sqrt(wsum / sum(wght * dhat^2))
  cini <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    cini[i, j] <- cini[j, i] <- dhat[k]
  }
  rc <- apply(cini, 1, mean)
  rm <- mean(cini)
  cini <- -(cini - outer(rc, rc, "+") + rm) /  2
  return(cini)
}

smacofBound <- function(theData, quick = FALSE) {
  ndat <- theData$ndat
  nobj <- theData$nobj
  iind <- theData$iind
  jind <- theData$jind
  wght <- theData$weights
  if (quick == 1) {
    return(2 * nobj * max(wght))
  }
  if (quick == 2) {
    return(4 * sum(wght))
  }
  asum <- 4 * diag(ndat)
  for (k in 1:(ndat - 1)) {
    kind <- c(iind[k], jind[k])
    for (l in (k + 1):ndat) {
      lind <- c(iind[l], jind[l])
      asum[k, l] <- asum[l, k] <- sum(outer(kind, lind, "=="))
    }
  }
  asum <- asum * sqrt(outer(wght, wght))
  return(eigs_sym(asum, 1)$values)
}