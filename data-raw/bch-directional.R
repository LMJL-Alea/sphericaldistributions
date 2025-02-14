library(sphericaldistributions)
mu <- c(1, 0, 0)
kappa <- 50
n <- 10000000
microbenchmark::microbenchmark(
  spl0 = Directional::rvmf(n, mu, kappa),
  spl1 = as_vmf_sample(rvmf(n, mu, kappa)),
  times = 10L
)

n <- 100
spl <- as_vmf_sample(rvmf(n, mu, kappa))


system.time(
  a <- Directional::dvmf(spl, mu, kappa)
)
system.time(
  b <- dvmf(spl, mu, kappa)
)
