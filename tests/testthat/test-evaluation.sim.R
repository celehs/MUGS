testthat::test_that("evaluation.sim handles minimal valid input", {
  set.seed(123)

  U <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rownames(U) <- paste0("code_", 1:10)

  pairs.rel <- data.frame(
    row = c("code_1", "code_2", "code_3", "code_4", "code_5"),
    col = c("code_6", "code_7", "code_8", "code_9", "code_10")
  )

  result <- MUGS::evaluation.sim(pairs.rel, U)

  testthat::expect_type(result, "list")
  testthat::expect_true("n.rel" %in% names(result))
  testthat::expect_true("AUC.rel" %in% names(result))
  testthat::expect_true(is.numeric(result$n.rel))
  testthat::expect_true(is.numeric(result$AUC.rel))
  testthat::expect_true(result$AUC.rel >= 0 && result$AUC.rel <= 1)
})

testthat::test_that("evaluation.sim returns error for empty pairs.rel", {
  set.seed(123)

  U <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  rownames(U) <- paste0("code_", 1:100)

  pairs.rel <- data.frame(row = character(0), col = character(0))

  testthat::expect_error(MUGS::evaluation.sim(pairs.rel, U), "'response' must have two levels")
})

testthat::test_that("evaluation.sim handles another valid case", {
  set.seed(123)

  U <- matrix(rnorm(20), nrow = 5, ncol = 4)
  rownames(U) <- paste0("code_", 1:5)

  pairs.rel <- data.frame(
    row = c("code_1", "code_2"),
    col = c("code_3", "code_4")
  )

  result <- MUGS::evaluation.sim(pairs.rel, U)

  testthat::expect_type(result, "list")
  testthat::expect_true("n.rel" %in% names(result))
  testthat::expect_true("AUC.rel" %in% names(result))
  testthat::expect_true(is.numeric(result$n.rel))
  testthat::expect_true(is.numeric(result$AUC.rel))
  testthat::expect_true(result$AUC.rel >= 0 && result$AUC.rel <= 1)
})
