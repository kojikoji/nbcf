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

test_that("p.width.vec from calculatePwidthVec is correct", {
  expect_equal(sum(calculatePwidthVec(seq(0, 1, by=0.1))),
               1)
                                  
})

test_that("p.dict from calculatePdict is correct", {
  expect_equal(dim(rcpp_calculate_p_dict(c(0.1, 0.5, 0.9), c(0.7, 1.2))),
               c(3, 2))
  expect_equal(
    rcpp_calculate_p_dict(c(0.1, 0.5, 0.9), c(0.7, 1.2)),
    matrix(c(0, 1, 2,
             0, 1, 2),
           ncol=2))
  p.dict <- rcpp_calculate_p_dict(seq(0, 1000)/1000, c(0.5, 2))
  expect_equal(
    p.dict[,1][1],
    0)
  expect_equal(
    p.dict[,1][1001],
    1000)
  expect_equal(
    p.dict[,2][1],
    0)
  expect_equal(
    p.dict[,2][1001],
    1000)

})


test_that("get_approximate_index is correct",{
  p.vec <- seq(0, 1000)/1000
  mod.p.vec <- calculatePvec(1001)
  approx.idx <- get_approximate_index(mod.p.vec, p.vec)
  expect_equal(
    approx.idx[300] + 1,
    which.min(abs(p.vec - mod.p.vec[300]))
  )
})
          


test_that("lnb probability from calculateLnbValues is correct", {
  lnb.values <- calculateLnbValues(c(10, 20), c(0.1, 0.5), 30)
  expect_equal(
    lnb.values[1, 1], log(dnbinom(10, 30, 1 - 0.1)))
  expect_equal(
    lnb.values[2, 1], log(dnbinom(10, 30, 1 - 0.5)))    
})

test_that("setupLnbDict work well", {
  lnb.dict <- setupLnbDict(1000, mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  expect_equal(
    dim(lnb.dict@values),
    c(300, 500))
  expect_gte(
    min(exp(lnb.dict@values)),
    0)
  expect_lte(
    max(exp(lnb.dict@values)),
    1)
  expect_equal(
    length(lnb.dict@count.dict),
    1001)
  expect_equal(
    dim(lnb.dict@p.dict),
    c(300, 800))
  expect_equal(
    nrow(lnb.dict@p.dict),
    length(lnb.dict@lprior.values))
  expect_gte(
    min(exp(lnb.dict@lprior.values)),
    0)
})
