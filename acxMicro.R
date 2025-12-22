library(microbenchmark)

h <- microbenchmark(
acx(xini, ord = 0, itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = 2, itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = 3, itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 1), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 2), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 3), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 0, 3), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 2, 3), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 3, 3), itmax = 1000, eps = 1e-6, verbose = FALSE))

microbenchmark(
acx(xini, ord = 0, itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = c(0, 3), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx2(xini, ord = c(0,3), itmax = 1000, eps = 1e-6, verbose = FALSE),
acx(xini, ord = 3, itmax = 1000, eps = 1e-6, verbose = FALSE),
acx2(xini, ord = 3, itmax = 1000, eps = 1e-6, verbose = FALSE))