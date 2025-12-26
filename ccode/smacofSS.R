dyn.load("smacofSS.so")

source("smacofDataUtilities.R")
source("smacofPlots.R")
source("smacofTorgerson.R")

smacofSS <- function(theData,
                     ndim = 2,
                     xinit = NULL,
                     ties = 1,
                     iord = 2,
                     safe = 1,
                     itmax = 1000,
                     eps = 1e-6,
                     digits = 10,
                     width = 15,
                     verbose = TRUE,
                     weighted = FALSE,
                     ordinal = FALSE) {
  if (is.null(xinit)) {
    xinit <- smacofTorgerson(theData, 2)$conf
  }
  xold <- xinit
  nord <- length(iord)
  nobj <- theData$nobj
  ndat <- theData$ndat
  itel <- 1
  iind <- theData$iind
  jind <- theData$jind
  dhat <- theData$delta
  wght <- theData$weights
  if (!weighted) {
    wght <- rep(1, ndat)
  }
  dhat <- dhat / sqrt(sum(wght * dhat^2))
  blks <- theData$blocks
  edis <- rep(0, ndat)
  for (k in 1:ndat) {
    i <- iind[k]
    j <- jind[k]
    edis[k] <- sqrt(sum((xold[i, ] - xold[j, ])^2))
  }
  sdd <- sum(wght * edis^2)
  sde <- sum(wght * dhat * edis)
  lbd <- sde / sdd
  edis <- lbd * edis
  xold <- lbd * xold
  sold <- sum(wght * (dhat - edis)^2)
  snew <- 0.0
  xold <- as.vector(xold)
  xnew <- xold
  h <- .C(
    "smacofSSEngine",
    nobj = as.integer(nobj),
    ndim = as.integer(ndim),
    ndat = as.integer(ndat),
    nord = as.integer(nord),
    safe = as.integer(safe),
    itel = as.integer(itel),
    ties = as.integer(ties),
    itmax = as.integer(itmax),
    digits = as.integer(digits),
    width = as.integer(width),
    verbose = as.integer(verbose),
    ordinal = as.integer(ordinal),
    weighted = as.integer(weighted),
    sold = as.double(sold),
    snew = as.double(snew),
    eps = as.double(eps),
    iind = as.integer(iind - 1),
    jind = as.integer(jind - 1),
    iord = as.integer(iord),
    blks = as.integer(blks),
    wght = as.double(wght),
    edis = as.double(edis),
    dhat = as.double(dhat),
    xold = as.double(xold),
    xnew = as.double(xnew)
  )
  result <- list(
    delta = theData$delta,
    dhat = h$dhat,
    confdist = h$edis,
    conf = matrix(h$xnew, nobj, ndim),
    weightmat = h$wght,
    stress = h$snew,
    ndim = ndim,
    init = xinit,
    niter = h$itel,
    nobj = nobj,
    iind = h$iind,
    jind = h$jind,
    weighted = weighted,
    ordinal = ordinal,
    ties = h$ties
  )
  class(result) <- c("smacofSSResult", "smacofSSUOResult")
  return(result)
}
