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

test_that("count.dict from calculateCountDict is correct", {
  expect_equal(calculateCountDict(as.integer(c(0, 3, 6))),
              c(0, 0, 3, 3, 3, 6, 6))
})
  
test_that("p.vec from calculatePvec is correct", {
  expect_equal(length(calculatePvec(10)),
               10)
  expect_equal(calculatePvec(10)[1],
               0)
  expect_equal(calculatePvec(10)[10],
               1)
})
