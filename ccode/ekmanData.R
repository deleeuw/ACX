data(ekman, package = "smacof")
ekman <- 1 - ekman
ekmanData <- makeMDSData(ekman)
ekmanLabels <- as.character(attr(ekman, "Labels"))

