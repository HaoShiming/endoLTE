# 创建测试文件
library(testthat)
library(endoLTE)

test_that("fbm function works correctly", {
  # Test basic functionality
  result <- fbm(H = 0.7, n = 100, T = 1, seed = 123)
  expect_type(result, "double")
  expect_length(result, 101)  # n + 1
  
  # Test error handling
  expect_error(fbm(H = 1.5, n = 100, T = 1))
  expect_error(fbm(H = 0.5, n = 0, T = 1))
})

test_that("generate_endoLTE_data produces correct structure", {
  data <- generate_endoLTE_data(N = 10, T = 20, seed = 123)
  
  expect_s3_class(data, "data.frame")
  expect_equal(nrow(data), 10 * 20)
  expect_true(all(c("unit_id", "time", "D", "y", "w", "z") %in% names(data)))
})

test_that("estimate_single_series returns valid output", {
  # Create simple test data
  set.seed(123)
  n <- 50
  time <- 1:n
  D <- c(rep(0, 25), rep(1, 25))
  w <- rnorm(n)
  z <- rnorm(n)
  y <- 1 + 0.5 * w + 0.2 * z + 2 * D + rnorm(n)
  
  result <- estimate_single_series(y = y, D = D, z = z, w = w, time = time, t_post = 1)
  
  expect_type(result, "double")
  expect_length(result, 1)
})