library(testthat)
library(dplyr)
library(pROC)

# source("R/evaluation_sim.R") # Adjust the path as needed.

test_that("evaluation.sim handles minimal valid input", {
  set.seed(123)

  # Increase the number of codes and dimensions
  U <- matrix(rnorm(50), nrow = 10, ncol = 5)  # 10 codes, embedding size 5
  rownames(U) <- paste0("code_", 1:10)

  # Provide multiple distinct pairs ensuring a decent pool of codes:
  pairs.rel <- data.frame(
    row = c("code_1", "code_2", "code_3", "code_4", "code_5"),
    col = c("code_6", "code_7", "code_8", "code_9", "code_10")
  )

  result <- evaluation.sim(pairs.rel, U)

  # Now we expect a list with 'n.rel' and 'AUC.rel'
  expect_type(result, "list")
  expect_true("n.rel" %in% names(result))
  expect_true("AUC.rel" %in% names(result))
  expect_true(is.numeric(result$n.rel))
  expect_true(is.numeric(result$AUC.rel))
  expect_true(result$AUC.rel >= 0 && result$AUC.rel <= 1)
})

test_that("evaluation.sim returns error for empty pairs.rel", {
  set.seed(123)

  U <- matrix(rnorm(100 * 10), nrow = 100, ncol = 10)
  rownames(U) <- paste0("code_", 1:100)

  pairs.rel <- data.frame(row = character(0), col = character(0))

  # Expect an error since auc requires two classes
  expect_error(evaluation.sim(pairs.rel, U), "'response' must have two levels")
})

test_that("evaluation.sim handles minimal valid input", {
  set.seed(123)

  # Create a small U matrix
  U <- matrix(rnorm(20), nrow = 5, ncol = 4)
  rownames(U) <- paste0("code_", 1:5)

  # Provide multiple pairs to ensure both positive and negative classes are created
  # This gives the function a better chance to produce two classes for AUC calculation
  pairs.rel <- data.frame(
    row = c("code_1", "code_2"),
    col = c("code_3", "code_4")
  )

  result <- evaluation.sim(pairs.rel, U)

  # Checks
  expect_type(result, "list")
  expect_true("n.rel" %in% names(result))
  expect_true("AUC.rel" %in% names(result))
  expect_true(is.numeric(result$n.rel))
  expect_true(is.numeric(result$AUC.rel))
  expect_true(result$AUC.rel >= 0 && result$AUC.rel <= 1)
})
