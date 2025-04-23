testthat::test_that("CodeSiteEff_l2_par works with valid input", {
  set.seed(123)
  n1 <- 5
  n2 <- 5
  p <- 5

  # Generate random matrices for S.1 and S.2
  S.1 <- matrix(rnorm(n1 * p), nrow = n1, ncol = p)
  S.2 <- matrix(rnorm(n2 * p), nrow = n2, ncol = p)

  # Assign shared column names
  colnames(S.1) <- paste0("code_", 1:p)
  colnames(S.2) <- paste0("code_", 1:p)

  common_codes <- intersect(colnames(S.1), colnames(S.2))
  testthat::expect_true(length(common_codes) > 0, info = "No overlapping column names found")

  U.1 <- matrix(rnorm(n1 * p), nrow = n1, ncol = p)
  U.2 <- matrix(rnorm(n2 * p), nrow = n2, ncol = p)
  V.1 <- U.1
  V.2 <- U.2

  delta.int <- matrix(0, nrow = n1 + n2, ncol = p)
  rownames(delta.int) <- c(paste0("code_", 1:n1), paste0("code_", (n1 + 1):(n1 + n2)))

  lambda.delta <- 0.1
  n.common <- length(common_codes)
  n.core <- 2

  result <- MUGS::CodeSiteEff_l2_par(
    S.1 = S.1, S.2 = S.2, n1 = n1, n2 = n2, U.1 = U.1, U.2 = U.2,
    V.1 = V.1, V.2 = V.2, delta.int = delta.int, lambda.delta = lambda.delta,
    p = p, common_codes = common_codes, n.common = n.common, n.core = n.core
  )

  testthat::expect_type(result, "list")
  testthat::expect_true(all(c("delta", "V.1.new", "V.2.new") %in% names(result)))
  testthat::expect_equal(dim(result$delta), c(n1 + n2, p))
  testthat::expect_equal(dim(result$V.1.new), c(n1, p))
  testthat::expect_equal(dim(result$V.2.new), c(n2, p))
})
