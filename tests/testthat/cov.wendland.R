
# Check output datatype
#----------------------------------------------------------------------------
test_that("Retains input object type", {
  theta <- c(0.5, 1, 0, 0, 0.3)

  distance_atomic <- 0.5
  distance_vector <- seq(0, 1, 0.05)
  distance_matrix <- as.matrix(dist(data.frame(x=runif(30), y = runif(30))))
  distance_spam <- spam::as.spam(distance_matrix)

  expect_is(cov.wendland(distance_vector, theta = theta), "numeric")
  expect_is(cov.wendland(distance_matrix, theta = theta), "matrix")
  expect_is(cov.wendland(distance_spam, theta = theta), "spam")
})


# Check behaviour below epsilon and above range
#----------------------------------------------------------------------------
test_that("Boundary behaviour correct", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  expect_equal(cov.wendland(h = c(.Machine$double.eps, 0.5), theta = theta),
               c(1.3, 0))})

# Integration settings - QNG is default and therefore tested previously
#----------------------------------------------------------------------------

test_that("QNG integration works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  cov.args <- list(numint.abstol = 1e-3, numint.reltol = 1e-3)
  expect_equal(cov.wendland(h = 0.4, theta = theta, cov.args = cov.args) > 0, TRUE)
})

test_that("QAG integration works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  cov.args <- list(numint.subintervals = 50, numint.qag_key = 1)
  expect_equal(cov.wendland(h = 0.4, theta = theta, cov.args = cov.args) > 0, TRUE)})

test_that("QAGS integration works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  cov.args <- list(numint.subintervals = 100, numint.qag_key = 1)
  covFun <- covarianceFactory(cov.wendland,
                              cov.args = list(numint.subintervals = 100))
  expect_equal(cov.wendland(h = 0.4, theta = theta, cov.args = cov.args) > 0, TRUE)})

# Interpolation settings
#----------------------------------------------------------------------------

test_that("Linear interpolation works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  cov.args <- list(interp.method = "linear", interp.num_support = 50)
  expect_equal(cov.wendland(h = 0.4, theta = theta, cov.args = cov.args) > 0, TRUE)})

test_that("Cubic spline interpolation works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  cov.args <- list(interp.method = "linear", interp.num_support = 50)
  expect_equal(cov.wendland(h = 0.4, theta = theta, cov.args = cov.args) > 0, TRUE)})

test_that("Polynomial interpolation works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  cov.args <- list(interp.method = "polynomial", interp.num_support = 50)
  expect_equal(cov.wendland(h = 0.4, theta = theta, cov.args = cov.args) > 0, TRUE)})

