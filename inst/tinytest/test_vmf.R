mu <- c(1, 0, 0)
n <- 10

# k = 0 -------------------------------------------------------------------

kappa <- 0

rvals <- rvmf(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dvmf(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pvmf(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_vmf_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)

# k = 10 ------------------------------------------------------------------

kappa <- 10

rvals <- rvmf(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dvmf(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pvmf(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_vmf_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)

# k = 100 -----------------------------------------------------------------

kappa <- 100

rvals <- rvmf(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dvmf(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pvmf(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_vmf_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)

# k = 1000 ----------------------------------------------------------------

kappa <- 1000

rvals <- rvmf(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dvmf(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pvmf(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_vmf_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)
