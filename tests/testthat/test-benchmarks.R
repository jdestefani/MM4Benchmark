context("Functional tests")
library(forecast) #library is allowed in tests but not in package

#Set up
X <- EuStockMarkets
splitting_point <- round(2*dim(X)[1]/3)
X_train <- X[1:splitting_point,]
n_row <- nrow(X)
n_col <- ncol(X)
h <- 5
expected_scalar_dimension <- 1

multivariate_results <- multivariate_M4benchmarks(X_train,h)

# test_that("M4_benchmarks", {
#   x <- EuStockMarkets[,1]
#   results <- M4_benchmarks(x,h)
#   expect_true(is.matrix(multivariate_results$Naive))
#   expect_equal(dim(multivariate_results$Naive),c(h,n_col))
# })

test_that("Naive - Forecast", {
  expect_true(is.matrix(multivariate_results$Naive))
  expect_equal(dim(multivariate_results$Naive),c(h,n_col))
})

test_that("Naive - Time", {
  expect_true(is.numeric(multivariate_results$TimeNaive))
  expect_equal(length(multivariate_results$TimeNaive),expected_scalar_dimension)
})

test_that("Naive seasonal - Forecast", {
  expect_true(is.matrix(multivariate_results$NaiveSeasonal))
  expect_equal(dim(multivariate_results$NaiveSeasonal),c(h,n_col))
})

test_that("Naive seasonal - Time", {
  expect_true(is.numeric(multivariate_results$TimeNaiveSeasonal))
  expect_equal(length(multivariate_results$TimeNaiveSeasonal),expected_scalar_dimension)
})

test_that("Naive 2 - Forecast", {
  expect_true(is.matrix(multivariate_results$Naive2))
  expect_equal(dim(multivariate_results$Naive2),c(h,n_col))
})

test_that("Naive 2 - Time", {
  expect_true(is.numeric(multivariate_results$TimeNaive2))
  expect_equal(length(multivariate_results$TimeNaive2),expected_scalar_dimension)
})

test_that("Exponential Smoothing - Forecast", {
  expect_true(is.matrix(multivariate_results$SimpleES))
  expect_equal(dim(multivariate_results$SimpleES),c(h,n_col))
})

test_that("Exponential Smoothing - Time", {
  expect_true(is.numeric(multivariate_results$TimeSimpleES))
  expect_equal(length(multivariate_results$TimeSimpleES),expected_scalar_dimension)
})

test_that("Holt Winters - Forecast", {
  expect_true(is.matrix(multivariate_results$HoltWinters))
  expect_equal(dim(multivariate_results$HoltWinters),c(h,n_col))
})

test_that("Holt Winters - Time", {
  expect_true(is.numeric(multivariate_results$TimeHoltWinters))
  expect_equal(length(multivariate_results$TimeHoltWinters),expected_scalar_dimension)
})

test_that("Holt Winters Damped - Forecast", {
  expect_true(is.matrix(multivariate_results$HoltWintersDamped))
  expect_equal(dim(multivariate_results$HoltWintersDamped),c(h,n_col))
})

test_that("Holt Winters Damped - Time", {
  expect_true(is.numeric(multivariate_results$TimeHoltWintersDamped))
  expect_equal(length(multivariate_results$TimeHoltWintersDamped),expected_scalar_dimension)
})

test_that("Theta - Forecast", {
  expect_true(is.matrix(multivariate_results$Theta))
  expect_equal(dim(multivariate_results$Theta),c(h,n_col))
})

test_that("Theta - Time", {
  expect_true(is.numeric(multivariate_results$TimeTheta))
  expect_equal(length(multivariate_results$TimeTheta),expected_scalar_dimension)
})

test_that("Combined - Forecast", {
  expect_true(is.matrix(multivariate_results$Comb))
  expect_equal(dim(multivariate_results$Comb),c(h,n_col))
})

test_that("Combined - Time", {
  expect_true(is.numeric(multivariate_results$TimeComb))
  expect_equal(length(multivariate_results$TimeComb),expected_scalar_dimension)
})
