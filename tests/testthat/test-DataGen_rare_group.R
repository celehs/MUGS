# test-DataGen_rare_group.R

testthat::test_that("DataGen_rare_group returns correct output structure and dimensions", {
  set.seed(123)

  # Parameters
  seed <- 1
  p <- 5
  n1 <- 100
  n2 <- 100
  n.common <- 50
  n.group <- 30
  sigma.eps.1 <- 1
  sigma.eps.2 <- 3
  ratio.delta <- 0.05
  network.k <- 5
  rho.beta <- 0.5
  rho.U0 <- 0.4
  rho.delta <- 0.7
  sigma.rare <- 10
  n.rare <- 20
  group.size <- 5

  # Run function
  output <- MUGS::DataGen_rare_group(
    seed = seed, p = p, n1 = n1, n2 = n2, n.common = n.common,
    n.group = n.group, sigma.eps.1 = sigma.eps.1, sigma.eps.2 = sigma.eps.2,
    ratio.delta = ratio.delta, network.k = network.k, rho.beta = rho.beta,
    rho.U0 = rho.U0, rho.delta = rho.delta, sigma.rare = sigma.rare,
    n.rare = n.rare, group.size = group.size
  )

  # Ensure output is a list
  expect_type(output, "list")

  # Expected component names
  expected_names <- c(
    "delta1", "delta2", "u.1", "u.2",
    "S.1", "S.2", "S.1.0", "S.2.0",
    "X.group.source", "X.group.target",
    "pairs.rel.CV", "pairs.rel.EV"
  )
  expect_true(all(expected_names %in% names(output)))

  # Dimensions
  expect_equal(dim(output$delta1), c(n1, p))
  expect_equal(dim(output$delta2), c(n2, p))
  expect_equal(dim(output$u.1), c(n1, p))
  expect_equal(dim(output$u.2), c(n2, p))
  expect_equal(dim(output$S.1), c(n1, n1))
  expect_equal(dim(output$S.2), c(n2, n2))
  expect_equal(dim(output$S.1.0), c(n1, n1))
  expect_equal(dim(output$S.2.0), c(n2, n2))
  expect_equal(dim(output$X.group.source), c(n1, n.group))
  expect_equal(dim(output$X.group.target), c(n2, n.group))

  # Data types
  expect_s3_class(output$pairs.rel.CV, "data.frame")
  expect_s3_class(output$pairs.rel.EV, "data.frame")
  expect_true(all(c("row", "col", "type") %in% colnames(output$pairs.rel.CV)))
  expect_true(all(c("row", "col", "type") %in% colnames(output$pairs.rel.EV)))

  # Check for missing values
  expect_false(anyNA(output$delta1))
  expect_false(anyNA(output$delta2))

  # Check symmetry
  expect_equal(output$S.1, t(output$S.1))
  expect_equal(output$S.2, t(output$S.2))
})
