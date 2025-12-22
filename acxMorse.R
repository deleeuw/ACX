data(morse, package = "smacof")
delta <- as.matrix(morse)
wmat <- 1 - diag(36)
vmat <- -wmat
diag(vmat) <- -rowSums(vmat)
vinv <- solve(vmat + 1/36) - 1/36

