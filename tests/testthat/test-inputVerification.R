context("Input Verification")

X <- matrix(AirPassengers,ncol = 12)
h <- 5
n_row <- dim(X)[1]
n_col <- dim(X)[2]
expected_length_scalar <- 1

test_that("input_check_multivariate - X not matrix", {
  X <- AirPassengers[1:10]
  h <- 5
  expect_error(input_check_multivariate(X,h),
               error_X_type_string,
               fixed=T)
})

test_that("input_check_multivariate - h not numeric", {
  X <- matrix(AirPassengers,ncol = 12)
  h <- "5"
  expect_error(input_check_multivariate(X,h),
               error_h_type_string,
               fixed=T)
})

test_that("input_check_univariate - X not matrix", {
  x <- c("1","2","3")
  h <- 5
  expect_error(input_check_univariate(x,h),
               error_x_type_string,
               fixed=T)
})

test_that("input_check_univariate - h not numeric", {
  x <- AirPassengers[1:10]
  h <- "5"
  expect_error(input_check_univariate(x,h),
               error_h_type_string,
               fixed=T)
})

test_that("input_check_level - level not numeric", {
  level <- c("1","2","3")
  expect_error(input_check_level(level),
               error_level_type_string,
               fixed=T)
})

test_that("input_check_level - level wrong length", {
  level <- c(1)
  expect_error(input_check_level(level),
               error_level_length_string,
               fixed=T)
})

test_that("input_check_level - level not numeric", {
  level <- c(90,85)
  expect_error(input_check_level(level),
               error_level_order_string,
               fixed=T)
})
