context("test-setupLnbDict")

test_that("length of count.vec from calculateCountVec is correct", {
  
  expect_equal(length(calculateCountVec(1000, 500, continuous.upper = 300)),
               500)
  expect_equal(length(calculateCountVec(100, 500, continuous.upper = 300)),
               101)
  expect_equal(length(calculateCountVec(1000, 200, continuous.upper = 300)),
               200)
})

test_that("count.vec from calculateCountVec is correct", {
  expect_equal(calculateCountVec(8, 5, 0), c(0, 1, 2, 4, 8))
  expect_equal(calculateCountVec(9, 4+1, 2), c(0, 1, 2, 3, 9))
})
  
