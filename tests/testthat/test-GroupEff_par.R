library(testthat)
library(parallel)
library(doSNOW)
library(glmnet)
library(Matrix)

testthat::test_that("GroupEff_par works with multiple groups in both sites", {
  set.seed(123)

  # Parameters
  n.MGB <- 10
  n.BCH <- 10
  p <- 5
  n.group <- 4
  n.core <- 2
  lambda <- 0

  # Create name list
  name.list <- paste0("code_", 1:(n.MGB + n.BCH))

  # MGB and BCH codes
  MGB.codes <- name.list[1:n.MGB]
  BCH.codes <- name.list[(n.MGB + 1):(n.MGB + n.BCH)]

  # S.MGB and S.BCH
  S.MGB <- matrix(rnorm(n.MGB^2), nrow = n.MGB)
  S.MGB <- (S.MGB + t(S.MGB)) / 2
  rownames(S.MGB) <- MGB.codes
  colnames(S.MGB) <- MGB.codes

  S.BCH <- matrix(rnorm(n.BCH^2), nrow = n.BCH)
  S.BCH <- (S.BCH + t(S.BCH)) / 2
  rownames(S.BCH) <- BCH.codes
  colnames(S.BCH) <- BCH.codes

  # U and V
  U.MGB <- matrix(rnorm(n.MGB * p), nrow = n.MGB)
  U.BCH <- matrix(rnorm(n.BCH * p), nrow = n.BCH)
  V.MGB <- matrix(rnorm(n.MGB * p), nrow = n.MGB)
  V.BCH <- matrix(rnorm(n.BCH * p), nrow = n.BCH)

  # Define groups
  group.names <- c("group_1", "group_2", "group_3", "group_4")

  X.MGB.group <- matrix(0, nrow = n.MGB, ncol = n.group,
                        dimnames = list(MGB.codes, group.names))
  X.MGB.group[, "group_1"] <- 1
  X.MGB.group[n.MGB - 1, "group_3"] <- 1
  X.MGB.group[n.MGB, "group_4"] <- 1

  X.BCH.group <- matrix(0, nrow = n.BCH, ncol = n.group,
                        dimnames = list(BCH.codes, group.names))
  X.BCH.group[, "group_2"] <- 1
  X.BCH.group[n.BCH - 1, "group_3"] <- 1
  X.BCH.group[n.BCH, "group_4"] <- 1

  beta.int <- matrix(0, nrow = n.group, ncol = p,
                     dimnames = list(group.names, NULL))

  # Run the function
  result <- MUGS::GroupEff_par(
    S.MGB = S.MGB, S.BCH = S.BCH, n.MGB = n.MGB, n.BCH = n.BCH,
    U.MGB = U.MGB, U.BCH = U.BCH, V.MGB = V.MGB, V.BCH = V.BCH,
    X.MGB.group = X.MGB.group, X.BCH.group = X.BCH.group,
    n.group = n.group, name.list = name.list,
    beta.int = beta.int, lambda = lambda, p = p, n.core = n.core
  )

  # Checks
  testthat::expect_type(result, "list")
  testthat::expect_true(all(c("beta", "dif_F", "V.MGB.new", "V.BCH.new") %in% names(result)))
  testthat::expect_equal(dim(result$beta), c(n.group, p))
  testthat::expect_true(is.numeric(result$dif_F))
  testthat::expect_equal(dim(result$V.MGB.new), c(n.MGB, p))
  testthat::expect_equal(dim(result$V.BCH.new), c(n.BCH, p))
})
