n <- 10

rvals <- runisph(n)
expect_equal(nrow(rvals), n)
expect_equal(ncol(rvals), 3)

dvals <- dunisph(rvals)
expect_equal(length(dvals), n)

pvals <- punisph(rvals)
expect_equal(length(pvals), n)
