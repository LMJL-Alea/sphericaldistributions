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

mval <- mean(as_watson_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)

# k = 10 ------------------------------------------------------------------

kappa <- 10

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_watson_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)

# k = 100 -----------------------------------------------------------------

kappa <- 100

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_watson_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)

# k = 1000 ----------------------------------------------------------------

kappa <- 1000

rvals <- rwatson(n, mu, kappa)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dwatson(rvals, mu, kappa)
expect_equal(length(dvals), n)

pvals <- pwatson(rvals, mu, kappa)
expect_equal(length(pvals), n)

mval <- mean(as_watson_sample(rvals))
expect_equal(nrow(mval), 1)
expect_equal(ncol(mval), 3)
