
# Check output datatype
#----------------------------------------------------------------------------
test_that("Retains input object type", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  distance_vector <- seq(0, 1, 0.05)
  distance_matrix <- as.matrix(dist(data.frame(x=runif(30), y = runif(30))))
  distance_spam <- as.spam(distance_matrix)
  covFun <- covarianceFactory(cov.wendland)
  expect_type(covFun(0.5, theta = theta), "double")
  expect_type(covFun(distance_vector, theta = theta), "double")
  expect_type(covFun(distance_matrix, theta = theta), "double")
  expect_s4_class(covFun(distance_spam, theta = theta), "spam")
})


# Check behaviour below epsilon and above range
#----------------------------------------------------------------------------
test_that("Boundary behaviour", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland)
  expect_equal(covFun(h = .Machine$double.eps, theta = theta), 1.3)
  expect_equal(covFun(0.5, theta = theta), 0)})


# Fixed parameter tests
#----------------------------------------------------------------------------
test_that("Fixed range functionality works", {
  theta <- c(1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland,
                              cov.args = list(fixed_range_value = 0.7))
  expect_equal(covFun(h = 0.8, theta = theta), 0)})

test_that("Fixed nugget functionality works", {
  theta <- c(0.5, 1, 0, 0)
  covFun <- covarianceFactory(cov.wendland,
                              cov.args = list(fixed_nugget_value = 0))
  expect_equal(covFun(h = 0, theta = theta), 1)})


# Integration settings - QNG is default and therefore tested previously
#----------------------------------------------------------------------------

test_that("QAGS integration works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland,
                              cov.args = list(numint.subintervals = 100))
  expect_equal(covFun(h = 0.4, theta = theta) > 0, TRUE)})

test_that("QAG integration works", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland,
                              cov.args = list(numint.subintervals = 100,
                                              numint.qag_key = 1))
  expect_equal(covFun(h = 0.4, theta = theta) > 0, TRUE)})


# Interpolation settings
#----------------------------------------------------------------------------

test_that("Linear interpolation methods work", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland,
                              list(interp.method = "linear",
                                   interp.num_support = 50))
  expect_equal(covFun(h = 0.4, theta = theta) > 0, TRUE)})

test_that("Cubic spline interpolation methods work", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland,
                              list(interp.method = "cspline",
                                   interp.num_support = 50))
  expect_equal(covFun(h = 0.4, theta = theta) > 0, TRUE)})

test_that("Polynomial interpolation methods work", {
  theta <- c(0.5, 1, 0, 0, 0.3)
  covFun <- covarianceFactory(cov.wendland,
                              list(interp.method = "polynomial",
                                   interp.num_support = 50))
  expect_equal(covFun(h = 0.4, theta = theta) > 0, TRUE)})

