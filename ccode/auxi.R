auxi <- function(dat, ini) {
  ndat <- length(dat$delta)
  nobj <- nrow(ini)
  ndim <- ncol(ini)
  iind <- dat$iind
  jind <- dat$jind
  xupd <- matrix(0, nobj, ndim)
  bmat <- matrix(0, nobj, nobj)
  for (k in 1:ndat) {
    is <- iind[k]
    js <- jind[k]
    elem <- dat$delta[k] / sqrt(sum((ini[is, ] - ini[js, ])^2)) 
    bmat[is, js] <- bmat[js, is] <- -elem
    for (s in 1:ndim) {
      add <- elem * (ini[is, s] - ini[js, s])
      xupd[is, s] <- xupd[is, s] + add
      xupd[js, s] <- xupd[js, s] - add
    }
  }
  diag(bmat) <- -rowSums(bmat)
 return(list(xupd / nobj, bmat %*% ini / nobj))
}
