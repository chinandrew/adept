
context("Testing pmax_max_cpp.cpp functions")

test_that("pmax_max_cpp is corrrect.", {
  x = 1:10
  y = c(5,0,5,0,5,0,5,0,5,0)
  z = c(12,0,0,0,0,0,0,0,0,12)
  inputs = list(x,y,z)
  output = pmax_max_cpp(inputs)
  expect_equal(output$pmax,
               c(12,2,5,4,5,6,7,8,9,12))
  expect_equal(output$idx,
               c(3,1,2,1,1,1,1,1,1,3))
})


test_that("pmax_max_cpp is agrees with previous implementation.", {
  x = 1:10
  y = c(5,0,5,0,5,0,5,0,5,0)
  z = c(12,0,0,0,0,0,0,0,0,12)
  inputs = list(x,y,z)
  output = pmax_max_cpp(inputs)
  expect_equal(output$pmax,
               do.call(pmax, inputs))
  expect_equal(output$idx,
               max.col(t(do.call(rbind, inputs)), ties.method = "first"))
})
