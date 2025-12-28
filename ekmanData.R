source("smacofDataUtilities.R")

data(ekman, package = "smacof")
ekmanDist <- 1 - ekman
ekmanMatrix <- as.matrix(ekmanDist)
ekmanData <- makeMDSData(ekmanDist)
ekmanLabels <- as.character(attr(ekman, "Labels"))

