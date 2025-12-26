source("smacofSS.R")

x <- matrix(c(0,0,1,0,1,1,0,1),4,2,byrow = TRUE)
small <- dist(x)
smallData <- makeMDSData(small)
small <- as.matrix(small)
small <- small / sqrt(sum(small^2))
x0 <- matrix(c(1:4,1,1,0,0),4,2)
x0 <- apply(x0, 2, function(x) x - mean(x))
d <- as.matrix(dist(x0))
lbd <- sum(small * d) / sum(d^2)
x0 <- x0 * lbd

d <- as.matrix(dist(x0))
print(sum((small - d)^2))
b <- -small / (d + diag(4))
diag(b)<- -rowSums(b)
x1 <- (b / 4) %*% x0

d <- as.matrix(dist(x1))
print(sum((small - d)^2))
b <- -small / (d + diag(4))
diag(b)<- -rowSums(b)
x2 <- (b / 4) %*% x1

d <- as.matrix(dist(x2))
print(sum((small - d)^2))
b <- -small / (d + diag(4))
diag(b)<- -rowSums(b)
x3 <- (b / 4) %*% x2

d0 <- x0
d1 <- x1 - x0
d2 <- x2 - 2 * x1 + x0
d3 <- x3 - 3 * x2 + 3 * x1 - x0

s1 <- abs(sum(d0 * d1) / sum(d1^2))
s2 <- abs(sum(d1 * d2) / sum(d2^2))
s3 <- abs(sum(d2 * d3) / sum(d3^2))

y1 <- d0 + s1 * d1
d <- as.matrix(dist(y1))
print(sum((small - d)^2))
y2 <- d0 + 2 * s2 * d1 + (s2 ^ 2) * d2
d <- as.matrix(dist(y2))
print(sum((small - d)^2))
y3 <- d0 + 3 * s3 * d1 + 3 * (s3 ^ 2) * d2 + (s3 ^ 3) * d3
d <- as.matrix(dist(y3))
print(sum((small - d)^2))
