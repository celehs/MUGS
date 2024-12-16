library(testthat)
library(MASS)
library(mvtnorm)
library(utils)
library(stats)
library(fastDummies)
library(Matrix)

test_that("DataGen_rare_group returns the correct output structure with updated parameters", {
  set.seed(123)

  # testing parameters
  seed <- 1
  p <- 100
  n1 <- 2000
  n2 <- 2000
  n.common <- 1000
  n.group <- 400
  sigma.eps.1 <- 5
  sigma.eps.2 <- 20
  sigma.rare <- 80
  ratio.delta <- 0.05
  network.k <- 50
  rho.beta <- 0.4
  rho.U0 <- 0.4
  rho.delta <- 0.95
  n.rare <- 700
  group.size <- 5

  # Run function
  output <- DataGen_rare_group(
    seed = seed, p = p, n1 = n1, n2 = n2, n.common = n.common, n.group = n.group,
    sigma.eps.1 = sigma.eps.1, sigma.eps.2 = sigma.eps.2, ratio.delta = ratio.delta,
    network.k = network.k, rho.beta = rho.beta, rho.U0 = rho.U0, rho.delta = rho.delta,
    sigma.rare = sigma.rare, n.rare = n.rare, group.size = group.size
  )

  # Check that output is a list
  expect_type(output, "list")

  # Check required list elements
  expected_names <- c("delta1", "delta2", "u.1", "u.2", "S.1", "S.2", "S.1.0", "S.2.0",
                      "X.group.source", "X.group.target", "pairs.rel.CV", "pairs.rel.EV")
  expect_true(all(expected_names %in% names(output)))

  # Check dimensions based on updated parameters:
  # delta1: n1 x p
  expect_equal(dim(output$delta1), c(n1, p))
  # delta2: n2 x p
  expect_equal(dim(output$delta2), c(n2, p))

  # u.1: n1 x p
  expect_equal(dim(output$u.1), c(n1, p))
  # u.2: n2 x p
  expect_equal(dim(output$u.2), c(n2, p))

  # S.1: n1 x n1
  expect_equal(dim(output$S.1), c(n1, n1))
  # S.2: n2 x n2
  expect_equal(dim(output$S.2), c(n2, n2))

  # S.1.0: n1 x n1
  expect_equal(dim(output$S.1.0), c(n1, n1))
  # S.2.0: n2 x n2
  expect_equal(dim(output$S.2.0), c(n2, n2))

  # X.group.source: n1 x n.group
  expect_equal(dim(output$X.group.source), c(n1, n.group))
  # X.group.target: (n2 - n.common) x n.group
  expect_equal(dim(output$X.group.target), c(2000, n.group))

  # pairs.rel.CV and pairs.rel.EV should be data frames
  expect_s3_class(output$pairs.rel.CV, "data.frame")
  expect_s3_class(output$pairs.rel.EV, "data.frame")

  # Check that pairs.rel.CV and pairs.rel.EV have columns "row", "col", and "type"
  expect_true(all(c("row", "col", "type") %in% colnames(output$pairs.rel.CV)))
  expect_true(all(c("row", "col", "type") %in% colnames(output$pairs.rel.EV)))

  # Basic sanity checks:
  # Ensure no NA values in delta matrices
  expect_false(anyNA(output$delta1))
  expect_false(anyNA(output$delta2))

  # Ensure S.1 and S.2 are symmetric
  expect_equal(output$S.1, t(output$S.1))
  expect_equal(output$S.2, t(output$S.2))
})
