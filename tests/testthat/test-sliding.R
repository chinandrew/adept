context("Testing sliding cov and cor functions.")

test_that("sliding_cov_fast agrees with dvmisc::sliding_cov", {
  long <- c(1, 5, 2, 7, 4, 5)
  short <- c(-1, 0, 1)
  expect_equal(dvmisc::sliding_cov(short, long),
               sliding_cov_fast(short, long))
})
