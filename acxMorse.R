source("smacofDataUtilities.R")

data(morse, package = "smacof")
morseDist <- morse
morseMatrix <- as.matrix(morse)
morseData <- makeMDSData(morseDist)

