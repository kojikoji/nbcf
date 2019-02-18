context("test calculateMpa.R")

test_that("traceGetShatl work well",{
  shatl <- c(1, 3, 5, 5, 5)
  taul <- traceGetShatl(shatl, c(0))
  expect_equal(
    taul,
    c(1, 3)
  )
})

test_that("calculateMap work well",{
  lpst <- matrix(runif(100, -10, -1), nrow=10)
  lambda <- 0.01
   map.change.point <- calculateMap(lpst, lambda)
  expect_lte(
    length(map.change.point),
    10)
})
