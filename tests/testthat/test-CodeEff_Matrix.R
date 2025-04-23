testthat::test_that("CodeEff_Matrix works as expected", {
  set.seed(1)

  n1 <- 100
  n2 <- 80
  p <- 10

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

  common_codes <- intersect(rownames(U.1), rownames(U.2))
  zeta.int <- matrix(0, n1 + n2, p)
  rownames(zeta.int) <- c(rownames(U.1), rownames(U.2))

  lambda <- 1

  output <- MUGS::CodeEff_Matrix(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda, p)

  testthat::expect_type(output, "list")
  testthat::expect_named(output, c("zeta", "dif_F", "V.1.new", "V.2.new"))

  testthat::expect_equal(dim(output$zeta), c(n1 + n2, p))
  testthat::expect_equal(dim(output$V.1.new), c(n1, p))
  testthat::expect_equal(dim(output$V.2.new), c(n2, p))

  testthat::expect_true(is.numeric(output$dif_F))
  testthat::expect_true(output$dif_F >= 0)

  testthat::expect_true(all(is.finite(output$zeta)))
  testthat::expect_true(all(is.finite(output$V.1.new)))
  testthat::expect_true(all(is.finite(output$V.2.new)))
})

testthat::test_that("CodeEff_Matrix handles no overlap in common codes", {
  set.seed(123)

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

  common_codes <- character(0)

  output <- MUGS::CodeEff_Matrix(S.1, S.2, n1, n2, U.1, U.2, V.1, V.2, common_codes, zeta.int, lambda = 1, p)

  testthat::expect_equal(dim(output$zeta), c(n1 + n2, p))
  testthat::expect_equal(dim(output$V.1.new), c(n1, p))
  testthat::expect_equal(dim(output$V.2.new), c(n2, p))
  testthat::expect_true(is.numeric(output$dif_F))
})
