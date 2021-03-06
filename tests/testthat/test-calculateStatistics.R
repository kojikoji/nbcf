context("test-setupLpst")


test_that("log_sum_exp works", {
  logSumExp <- function(x) log(sum(exp(x - max(x)))) + max(x)
  vec <- c(-2, -12, -8)
  expect_equal(
    log_sum_exp(vec, c(1, 1, 1)),
    logSumExp(vec))
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

test_that("rcpp_calculate_lpst_g work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  t.grids <- seq(0, 10)*80
  lpst_g <- rcpp_calculate_lpst_g(count.mat[1, ], lnb.dict@values, lnb.dict@lprior.values,
                 lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
  expect_equal(sum(is.nan(lpst_g)),
              0)
  expect_equal(dim(lpst_g),
              c(10, 10))
})

test_that("calculate_lpst_list work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  t.grids <- seq(0, 10)*80
  lpst_list <- rcpp_calculate_lpst_list(count.mat, lnb.dict@values, lnb.dict@lprior.values,
                 lnb.dict@count.dict, lnb.dict@p.dict, lnb.dict@p.width.vec, t.grids)
  expect_equal(length(lpst_list),
               6)
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

test_that("calculateLpstVar work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  t.grids <- seq(0, 10)*80
  mean.bias.vec <- rep(1, 800)
  rownames(count.mat) <- as.character(seq(nrow(count.mat)))
  lpst <- calculateVarLpstList(count.mat, mean.bias.vec, t.grids)
  expect_equal(sum(is.nan(lpst[[1]])),
              0)
  expect_equal(dim(lpst[[1]]),
              c(10, 10))
})


test_that("calculateLpst work well", {
  count.mat <- Matrix(rnbinom(6*800, 30, 0.995), nrow=6)
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec=runif(800, 0.8, 1.2),
                           r=30, alpha=2, beta=2, count.res=500, p.res=300)
  t.grids <- seq(0, 10)*80
  one.var.lpst <- calculateOneVariateLpst(lnb.dict, count.mat[1, ], t.grids)
  expect_equal(dim(one.var.lpst),
              c(10, 10))
})


test_that("calculateQt work well", {
  lpst <- matrix(runif(10*10, -1, -0.01), nrow=10)
  lambda <- 0.1
  lqt <- calculateLqt(lpst, lambda)
  expect_equal(
    length(lqt),
    10)
})
