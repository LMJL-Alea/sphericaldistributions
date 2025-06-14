mu <- c(1, 0, 0)
n <- 10

# k = 0 -------------------------------------------------------------------

kappa <- 0

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(rvals)
expect_length(mval, 3L)

# k = 10 ------------------------------------------------------------------

kappa <- 10

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(rvals)
expect_length(mval, 3L)

# k = 100 -----------------------------------------------------------------

kappa <- 100

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(rvals)
expect_length(mval, 3L)

# k = 1000 ----------------------------------------------------------------

kappa <- 1000

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(rvals)
expect_length(mval, 3L)
