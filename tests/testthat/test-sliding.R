context("Testing sliding cov and cor functions.")

test_that("sliding_cov2 agrees with dvmisc::sliding_cov", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(1.5, 1, 2.1)
  expect_equal(dvmisc::sliding_cov(short, long),
               sliding_cov_fast(short, long))
})

test_that("sliding_cor2 agrees with dvmisc::sliding_cov", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(1.5, 1, 2.1)
  expect_equal(dvmisc::sliding_cor(short, long),
               sliding_cor_fast(short, long))
})
