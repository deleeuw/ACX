library(RSpectra)

torgerson <- function(delta, ndim = 2) {
  dd <- delta^2
  rd <- apply(dd, 1, mean)
  md <- mean(dd)
  ed <- -(dd - outer(rd, rd, "+") + md) / 2
  ev <- eigs_sym(ed, ndim)
  return(ev$vectors[, 1:ndim] %*% diag(sqrt(ev$values[1:ndim])))
}

guttman <- function(x, delta, wght, vinv) {
  n <- nrow(x)
  dmat <- as.matrix(dist(x))
  emat <- -wght * delta / (dmat + diag(n))
  diag(emat) <- -rowSums(emat)
  return(vinv %*% emat %*% x)
}

loss <- function(x, delta, wght) {
  dmat <- as.matrix(dist(x))
  return(sum(wght * (delta - dmat)^2) / 2)
}

mPrint <- function(x,
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


smacofACX <- function(delta,
                      wght = 1 - diag(nrow(delta)),
                      ndim = 2,
                      xinit = torgerson(delta, ndim),
                      strategy = 2,
                      width = 12,
                      digits = 10,
                      itmax = 1000,
                      eps = 1e-6,
                      safe = TRUE,
                      verbose = TRUE) {
  itel <- 1
  nobj <- nrow(delta)
  p <- length(strategy)
  vmat <- -wght
  diag(vmat) <- -rowSums(vmat)
  vinv <- solve(vmat + (1 / nobj)) - (1 / nobj)
  xold <- xinit
  fold <- loss(xold, delta, wght)
  repeat {
    pk <- strategy[(itel - 1) %% p + 1]
    d0 <- xold
    x1 <- guttman(xold, delta, wght, vinv)
    f1 <- loss(x1, delta, wght)
    if (pk == 0) {
      sigd <- 0
      xnew <- x1
    } else {
      d1 <- x1 - xold
      if (pk == 1) {
        sigd <- abs(sum(d1 * d0) / sum(d1^2))
        xnew <- d0 + sigd * d1
      } else {
        x2 <- guttman(x1, delta, wght, vinv)
        d2 <- x2 - 2 * x1 + xold
        sigd <- abs(sum(d2 * d1) / sum(d2^2))
        if (pk == 2) {
          xnew <- d0 + 2 * sigd * d1 + sigd^2 * d2
        } else {
          x3 <- guttman(x2, delta, wght, vinv)
          d3 <- x3 - 3 * x2 + 3 * x1 - xold
          if (pk == 3) {
            sigd <- abs(sum(d3 * d2)) / sum(d3^2)
            xnew <- d0 + 3 * sigd * d1 + 3 * sigd^2 * d2 + sigd^3 * d3
          }
        }
      }
    }
    fnew <- loss(xnew, delta, wght)
    if (safe && (f1 < fnew)) {
      fnew <- f1
      xnew <- x1
      sigd <- 0
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
        "sigd",
        formatC(
          sigd,
          format = "f",
          width = width,
          digits = digits
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
