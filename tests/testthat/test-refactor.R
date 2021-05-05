context("Refactoring of benchmarks - Equality between joint and individual functions")

#Set up
X <- AirPassengers
splitting_point <- round(2*length(X)/3)
X_train <- X[1:splitting_point]
h <- 5

joint_results <- M4_benchmarks(X_train,h)

test_that("Naive", {
  individual_results <- naiveBenchmark(X_train,h)
  expect_equal(joint_results$Naive,individual_results$forecasts)
  cat(paste("\nNaive\n","\nTime grouped:",joint_results$TimeNaive,"\nTime individual:",individual_results$time))
})

test_that("Naive seasonal", {
  individual_results <- naiveSeasonalBenchmark(X_train,h)
  expect_equal(joint_results$NaiveSeasonal,individual_results$forecasts)
  cat(paste("\n\nNaive seasonal\n","\nTime grouped:",joint_results$TimeNaiveSeasonal,"\nTime individual:",individual_results$time))
})

test_that("Naive 2 seasonal", {
  individual_results <- naive2Benchmark(X_train,h)
  expect_equal(joint_results$Naive2,individual_results$forecasts)
  cat(paste("\n\nNaive 2 seasonal\n","\nTime grouped:",joint_results$TimeNaive2,"\nTime individual:",individual_results$time))
})

test_that("Exponential Smoothing", {
  individual_results <- ESBenchmark(X_train,h)
  expect_equal(joint_results$SimpleES,individual_results$forecasts)
  cat(paste("\n\nES\n","\nTime grouped:",joint_results$TimeSimpleES,"\nTime individual:",individual_results$time))
})

test_that("Holt Winters", {
  individual_results <- HoltWintersBenchmark(X_train,h)
  expect_equal(joint_results$HoltWinters,individual_results$forecasts)
  cat(paste("\n\nHolt Winters\n","\nTime grouped:",joint_results$TimeHoltWinters,"\nTime individual:",individual_results$time))
})

test_that("Holt Winters Damped", {
  individual_results <- HoltWintersDampedBenchmark(X_train,h)
  expect_equal(joint_results$HoltWintersDamped,individual_results$forecasts)
  cat(paste("\n\nHolt Winters Damped\n","\nTime grouped:",joint_results$TimeHoltWintersDamped,"\nTime individual:",individual_results$time))
})

test_that("Theta properly implemented", {
  individual_results <- thetaBenchmark(X_train,h)
  expect_equal(joint_results$Theta,individual_results$forecasts)
  cat(paste("\n\nTheta\n","\nTime grouped:",joint_results$TimeTheta,"\nTime individual:",individual_results$time))
})

test_that("Combined", {
  individual_results <- combinedBenchmark(X_train,h)
  expect_equal(joint_results$Comb,individual_results$forecasts)
  cat(paste("\n\nCombined\n","\nTime grouped:",joint_results$TimeComb,"\nTime individual:",individual_results$time))
})
