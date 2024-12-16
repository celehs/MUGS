test_that("CodeEff_Matrix works as expected", {
  set.seed(1)

  # Simulate input data
  n1 <- 100   # Number of source site codes
  n2 <- 80    # Number of target site codes
  p <- 10     # Embedding length

  S.1 <- matrix(rnorm(n1 * n1), n1, n1)  # Source site SPPMI
  S.2 <- matrix(rnorm(n2 * n2), n2, n2)  # Target site SPPMI
  U.1 <- matrix(rnorm(n1 * p), n1, p)    # Left embeddings (source)
  U.2 <- matrix(rnorm(n2 * p), n2, p)    # Left embeddings (target)
  V.1 <- matrix(rnorm(n1 * p), n1, p)    # Right embeddings (source)
  V.2 <- matrix(rnorm(n2 * p), n2, p)    # Right embeddings (target)

  rownames(U.1) <- paste0("code:S_", 1:n1)
  rownames(U.2) <- paste0("code:T_", 1:n2)
  rownames(S.1) <- colnames(S.1) <- rownames(U.1)
  rownames(S.2) <- colnames(S.2) <- rownames(U.2)

  common_codes <- intersect(rownames(U.1), rownames(U.2))  # Find overlap
  zeta.int <- matrix(0, n1 + n2, p)  # Initialize zeta
  rownames(zeta.int) <- c(rownames(U.1), rownames(U.2))

  # Run CodeEff_Matrix
  lambda <- 1
  output <- CodeEff_Matrix(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda, p)

  # Test output structure
  expect_type(output, "list")
  expect_named(output, c("zeta", "dif_F", "V.1.new", "V.2.new"))

  # Test dimensions of zeta
  expect_equal(dim(output$zeta), c(n1 + n2, p))
  expect_equal(dim(output$V.1.new), c(n1, p))
  expect_equal(dim(output$V.2.new), c(n2, p))

  # Test Frobenius norm difference
  expect_true(is.numeric(output$dif_F))
  expect_true(output$dif_F >= 0)  # Should always be non-negative

  # Check that zeta values are reasonable
  expect_true(all(is.finite(output$zeta)))

  # Check updated embeddings
  expect_true(all(is.finite(output$V.1.new)))
  expect_true(all(is.finite(output$V.2.new)))
})

test_that("CodeEff_Matrix handles edge cases", {
  set.seed(123)

  # Edge case: Empty overlap (no common codes)
  n1 <- 10
  n2 <- 10
  p <- 5

  S.1 <- matrix(rnorm(n1 * n1), n1, n1)
  S.2 <- matrix(rnorm(n2 * n2), n2, n2)
  U.1 <- matrix(rnorm(n1 * p), n1, p)
  U.2 <- matrix(rnorm(n2 * p), n2, p)
  V.1 <- matrix(rnorm(n1 * p), n1, p)
  V.2 <- matrix(rnorm(n2 * p), n2, p)

  rownames(U.1) <- paste0("code:S_", 1:n1)
  rownames(U.2) <- paste0("code:T_", 1:n2)
  rownames(S.1) <- colnames(S.1) <- rownames(U.1)
  rownames(S.2) <- colnames(S.2) <- rownames(U.2)

  zeta.int <- matrix(0, n1 + n2, p)
  rownames(zeta.int) <- c(rownames(U.1), rownames(U.2))

  common_codes <- character(0)  # No common codes

  output <- CodeEff_Matrix(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda = 1, p)

  expect_equal(dim(output$zeta), c(n1 + n2, p))
  expect_equal(dim(output$V.1.new), c(n1, p))
  expect_equal(dim(output$V.2.new), c(n2, p))
  expect_true(is.numeric(output$dif_F))
})
