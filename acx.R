library(RSpectra)

torgerson <- function(delta, ndim = 2) {
  dd <- delta^2
  rd <- apply(dd, 1, mean)
  md <- mean(dd)
  ed <- -(dd - outer(rd, rd, "+") + md) / 2
  ev <- eigs_sym(ed, ndim)
  return(ev$vectors[, 1:ndim] %*% diag(sqrt(ev$values[1:ndim])))
}

guttman <- function(x) {
  n <- nrow(x)
  dmat <- as.matrix(dist(x))
  emat <- -wmat * delta / (dmat + diag(n))
  diag(emat) <- -rowSums(emat)
  return(vinv %*% emat %*% x)
}

loss <- function(x) {
  dmat <- as.matrix(dist(x))
  return(sum(wmat * (delta - dmat)^2))
}

smacofACX <- function(delta,
                      wght = 1 - diag(nrow(delta)),
                      ndim = 2,
                      xini = torgerson(delta, ndim),
                      ord = c(0, 3),
                      width = 12,
                      digits = 10,
                      itmax = 100,
                      eps = 1e-6,
                      mono = TRUE,
                      verbose = TRUE) {
  itel <- 0
  nobj <- nrow(deltat())
  p <- length(ord)
  vmat <- -wght
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  xold <- xini
  fold <- loss(xold)
  repeat {
    pk <- ord[itel %% p + 1]
    d0 <- xold
    x1 <- guttman(xold)
    x2 <- guttman(x1)
    f1 <- loss(x1)
    d1 <- x1 - xold
    d2 <- x2 - 2 * x1 + xold
    if (pk == 0) {
      xnew <- x1
    }
    if (pk == 1) {
      sigd <- abs(sum(d1 * d0) / sum(d1^2))
      xnew <- x1 + sigd * d1
    }
    if (pk == 2) {
      sigd <- abs(sum(d2 * d1) / sum(d2^2))
      xnew <- d0 + sigd * d1 + sigd^2 * d2
    }
    if (pk == 3) {
      x3 <- guttman(x2)
      d3 <- x3 - 3 * x2 + 3 * x1 - xold
      sigd <- abs(sum(d3 * d2)) / sum(d3^2)
      xnew <- d0 + 3 * sigd * d1 + 3 * sigd^2 * d2 + sigd^3 * d3
    }
    fnew <- loss(xnew)
    if (mono && (f1 < fnew)) {
      fnew <- f1
      xnew <- x1
      pk <- 0
    }
    difx <- max(abs(xold - xnew))
    if (verbose) {
      cat(
        "itel",
        formatC(itel, format = "d", width = 2),
        formatC(
          pk,
          format = "d",
          digits = 1,
          width = 1
        ),
        "difx",
        formatC(
          difx,
          format = "f",
          width = width,
          digits = digits
        ),
        "fold",
        formatC(
          fold,
          format = "f",
          width = width,
          digits = digits
        ),
        "fnew",
        formatC(
          fnew,
          format = "f",
          width = width,
          digits = digits
        ),
        "\n"
      )
    }
    if ((itel == itmax) || (difx < eps)) {
      break
    }
    itel <- itel + 1
    xold <- xnew
    fold <- fnew
  }
  return(list(x = xnew, itel = itel))
}
