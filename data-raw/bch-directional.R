library(sphericaldistributions)
mu <- c(1, 0, 0)
kappa <- 50
n <- 10000000
microbenchmark::microbenchmark(
  spl0 = Directional::rvmf(n, mu, kappa),
  spl1 = as_vmf_sample(rvmf(n, mu, kappa)),
  times = 10L
)

n <- 1000000
spl <- as_vmf_sample(rvmf(n, mu, kappa))
microbenchmark::microbenchmark(
  a = mean(spl),
  b = Directional::vmf.mle(spl, fast = FALSE, tol = 1e-15),
  c = Directional::vmf.mle(spl, fast = TRUE, tol = 1e-15)
)


system.time(
  a <- Directional::dvmf(spl, mu, kappa)
)
system.time(
  b <- dvmf(spl, mu, kappa)
)
