context("test-setupLpst")

test_that("test accumulate_colwise", {
  expect_equal(
    accumulate_colwise(
      matrix(c(1, 2, 3, 4, 5, 6), nrow=2)),
    matrix(c(1, 2, 1+3, 2+4, 1+3+5, 2+4+6), nrow=2))
})

test_that("log_sum_exp works", {
  logSumExp <- function(x) log(sum(exp(x - max(x)))) + max(x)
  vec <- c(-2, -12, -8)
  expect_equal(
    log_sum_exp(vec, c(1, 1, 1)),
    logSumExp(vec))
})

test_that("log_sum_exp works", {
  logSumExp <- function(x) log(sum(exp(x - max(x)))) + max(x)
  vec1 <- c(-2, -12, -8)
  vec2 <- c(-1, -9, -2)
  expect_equal(
    colwise_log_sum_exp(matrix(c(vec1, vec2), ncol=2), rep(1, 3)),
    c(logSumExp(vec1), logSumExp(vec2)))
})

test_that("calculate_pt_grid_lp work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  t.grids <- seq(0, 10)*80
  pt_grid_lp <- calculate_pt_grid_lp(count.mat[1,], lnb.dict@values,
                       lnb.dict@count.dict, lnb.dict@p.dict, t.grids)
  expect_equal(sum(is.nan(pt_grid_lp)),
              0)
  expect_equal(dim(pt_grid_lp),
               c(300, 10))
})

test_that("calculate_lpst_g work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  t.grids <- seq(0, 10)*80
  lpst_g <- calculate_lpst_g(count.mat[1, ], lnb.dict@values, lnb.dict@lprior.values,
                 lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
  expect_equal(sum(is.nan(lpst_g)),
              0)
  expect_equal(dim(lpst_g),
              c(10, 10))
})

test_that("calculateLpst work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  t.grids <- seq(0, 10)*80
  lpst <- calculateLpst(lnb.dict, count.mat, t.grids)
  expect_equal(sum(is.nan(lpst)),
              0)
  expect_equal(dim(lpst),
              c(10, 10))
})

