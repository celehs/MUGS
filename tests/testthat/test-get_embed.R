testthat::test_that("get_embed returns correct embedding with default parameters", {
  set.seed(123)
  n <- 100
  p <- 3000
  U <- matrix(rnorm(n * p), n, p)
  V <- matrix(rnorm(n * p), n, p)
  d_values <- sort(runif(p, min = 1, max = 10), decreasing = TRUE)
  names_vec <- paste0("code_", 1:n)

  mysvd <- list(u = U, v = V, d = d_values, names = names_vec)

  embed <- MUGS::get_embed(mysvd)

  testthat::expect_true(is.matrix(embed))
  testthat::expect_equal(nrow(embed), n)
  testthat::expect_equal(rownames(embed), names_vec)

  norms <- apply(embed, 1, function(x) sqrt(sum(x^2)))
  testthat::expect_true(all(abs(norms - 1) < 1e-8))
})

testthat::test_that("get_embed returns correct embedding with different d", {
  set.seed(123)
  n <- 50
  p <- 500
  U <- matrix(rnorm(n * p), n, p)
  V <- matrix(rnorm(n * p), n, p)
  d_values <- sort(runif(p, min = 1, max = 10), decreasing = TRUE)
  names_vec <- paste0("elem_", 1:n)

  mysvd <- list(u = U, v = V, d = d_values, names = names_vec)

  d_custom <- 100
  embed <- MUGS::get_embed(mysvd, d = d_custom, normalize = FALSE)

  testthat::expect_equal(nrow(embed), n)
  testthat::expect_true(ncol(embed) <= d_custom)

  norms <- apply(embed, 1, function(x) sqrt(sum(x^2)))
  testthat::expect_true(any(abs(norms - 1) > 1e-8))
})

testthat::test_that("get_embed handles cases where d > length(id)", {
  set.seed(123)
  n <- 20
  p <- 50
  U <- matrix(rnorm(n * p), n, p)
  V <- matrix(rnorm(n * p), n, p)
  d_values <- sort(runif(p, min = 1, max = 10), decreasing = TRUE)
  names_vec <- paste0("x_", 1:n)

  mysvd <- list(u = U, v = V, d = d_values, names = names_vec)

  embed <- MUGS::get_embed(mysvd, d = p + 100, normalize = TRUE)

  testthat::expect_true(ncol(embed) <= p)
  norms <- apply(embed, 1, function(x) sqrt(sum(x^2)))
  testthat::expect_true(all(abs(norms - 1) < 1e-8))
})
