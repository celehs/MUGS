library(testthat)

test_that("get_embed returns correct embedding with default parameters", {
  # Simulate SVD results
  set.seed(123)
  n <- 100    # number of rows
  p <- 3000   # number of singular values/vectors
  U <- matrix(rnorm(n * p), n, p)
  V <- matrix(rnorm(n * p), n, p)
  d_values <- sort(runif(p, min = 1, max = 10), decreasing = TRUE)
  names_vec <- paste0("code_", 1:n)

  # Create a fake mysvd object
  mysvd <- list(u = U, v = V, d = d_values, names = names_vec)

  # Default d=2000, normalize=TRUE
  embed <- get_embed(mysvd)

  # Checks
  expect_true(is.matrix(embed))
  # The embedding dimension should be min(d, length(id)) rows = n, columns <= d
  # We know by default d=2000, but id is determined by sign conditions.
  # Check row count matches number of names
  expect_equal(nrow(embed), n)
  # Check that rownames match mysvd$names
  expect_equal(rownames(embed), names_vec)

  # Check normalization: each row should have norm ~ 1
  norms <- apply(embed, 1, function(x) sqrt(sum(x^2)))
  expect_true(all(abs(norms - 1) < 1e-8))
})

test_that("get_embed returns correct embedding with different d", {
  set.seed(123)
  n <- 50
  p <- 500
  U <- matrix(rnorm(n * p), n, p)
  V <- matrix(rnorm(n * p), n, p)
  d_values <- sort(runif(p, min = 1, max = 10), decreasing = TRUE)
  names_vec <- paste0("elem_", 1:n)

  mysvd <- list(u = U, v = V, d = d_values, names = names_vec)

  # Choose a d smaller than p
  d_custom <- 100
  embed <- get_embed(mysvd, d = d_custom, normalize = FALSE)

  # Checks
  expect_equal(nrow(embed), n)
  # The number of columns should be at most d_custom, but can be less if id is shorter.
  # Since id is defined by sign conditions, test that it's <= d_custom.
  expect_true(ncol(embed) <= d_custom)

  # Check no normalization: row norms won't be 1
  norms <- apply(embed, 1, function(x) sqrt(sum(x^2)))
  expect_true(any(abs(norms - 1) > 1e-8)) # At least some deviation
})

test_that("get_embed handles cases where d > length(id)", {
  set.seed(123)
  n <- 20
  p <- 50
  U <- matrix(rnorm(n * p), n, p)
  V <- matrix(rnorm(n * p), n, p)
  d_values <- sort(runif(p, min = 1, max = 10), decreasing = TRUE)
  names_vec <- paste0("x_", 1:n)

  mysvd <- list(u = U, v = V, d = d_values, names = names_vec)

  # Set d larger than p
  embed <- get_embed(mysvd, d = p + 100, normalize = TRUE)

  # The embedding should not exceed p or the length of the id vector.
  expect_true(ncol(embed) <= p)
  # Check normalization again
  norms <- apply(embed, 1, function(x) sqrt(sum(x^2)))
  expect_true(all(abs(norms - 1) < 1e-8))
})
