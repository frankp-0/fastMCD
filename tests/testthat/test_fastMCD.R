test_that("Test cov", {
    set.seed(5500)
    sigma <- matrix(runif(10^2) * 2 - 1, 10)
    sigma <- t(sigma) %*% sigma
    X <- MASS::mvrnorm(1e4, mu = rep(0, 10), Sigma = sigma)
    X[seq(1, 1e4, 5), ] <- MASS::mvrnorm(2000, mu = rep(3, 10), Sigma = sigma)
    res <- fastMCD(X)
    expect_equal(res$cov, sigma, tolerance = 0.2)
})

test_that("Test center", {
    set.seed(80124)
    sigma <- matrix(runif(10^2) * 2 - 1, 10)
    sigma <- t(sigma) %*% sigma
    X <- MASS::mvrnorm(1e4, mu <- runif(10), Sigma = sigma)
    X[seq(1, 1e4, 5), ] <- MASS::mvrnorm(2000, mu = rep(3, 10), Sigma = sigma)
    res <- fastMCD(X)
    expect_gte(cor(mu, res$center), 0.9)
})
