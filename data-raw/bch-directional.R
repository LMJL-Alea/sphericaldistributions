library(sphericaldistributions)
mu <- c(1, 0, 0)
kappa <- 50
n <- 10000000
microbenchmark::microbenchmark(
  spl0 = Directional::rvmf(n, mu, kappa),
  spl1 = rvmf(n, mu, kappa),
  times = 10L
)

n <- 10000000
spl <- rvmf(n, mu, kappa)
bm <- bench::mark(
  a = mean(spl),
  b = Directional::vmf.mle(spl, fast = FALSE, tol = 1e-15)$mu,
  c = Directional::vmf.mle(spl, fast = TRUE, tol = 1e-15)$mu
)
