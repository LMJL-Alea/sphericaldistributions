x <- rwatson(100, c(1, 0, 0), 10)
plot(x)
expect_equal(nrow(x), 100L)
