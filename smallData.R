source("smacofSS.R")

smallMatrix <- matrix(c(
  0,1,2,4,7,11,
  1,0,3,5,8,12,
  2,3,0,6,9,13,
  4,5,6,0,10,14,
  7,8,9,10,0,15,
  11,12,13,14,15,0), 6, 6)
smallMatrix <- smallMatrix * sqrt(2 / sum(smallMatrix^2))
smallDist <- as.dist(smallMat)
smallData <- makeMDSData(smallDist)
xinit <- matrix(c(1,-1,1,-1,1, -1, 0,0,-1,-1, 1, 1),6,2)
dinit <- as.matrix(dist(xinit))
labda <- sum(smallMat * dinit) / sum(dinit ^ 2)
xinit <- labda * xinit
