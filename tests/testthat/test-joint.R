context("Joint test")

test_that("perfectSimulaton correctly mine change points", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 100
  t.vec <- seq(tnum)
  true.change.point <- 40
  before.change.vec <- rnbinom(6*(true.change.point), 30, 0.995)
  after.change.vec <- rnbinom(6*(tnum - true.change.point), 30, 0.95)
  count.mat <- Matrix(c(before.change.vec, after.change.vec), nrow=6, sparse=T)
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  nbcf <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                          alpha, beta, r, lambda,
                          p.res, count.res, t.res)
  ## check distance from tru points
  mean.dist.from.40 <- mean(unlist(purrr::map(nbcf@sim.change.point.list,
                                              ~ min(abs(.x - 40)))))
  expect_lt(
    mean.dist.from.40,
    10)
  ## check change points num
  median.change.point.num <- median(unlist(purrr::map(nbcf@sim.change.point.list,
                                                      ~ length(.x))))
  expect_equal(
    median.change.point.num,
    1)
})


test_that("perfectSimulaton correctly mine two wchange points", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 700
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  nbcf <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                          alpha, beta, r, lambda,
                          p.res, count.res, t.res)
  ## check distance from true points
  mean.dist.from.400 <- mean(unlist(purrr::map(nbcf@sim.change.point.list,
                                  ~ min(abs(.x - 400)))))
  expect_lt(
    mean.dist.from.400,
    10)
  mean.dist.from.700 <- mean(unlist(purrr::map(nbcf@sim.change.point.list,
                                  ~ min(abs(.x - 700)))))
  expect_lt(
    mean.dist.from.700,
    10)
  ## check change points num
  median.change.point.num <- median(unlist(purrr::map(nbcf@sim.change.point.list,
                                                      ~ length(.x))))
  expect_equal(
    median.change.point.num,
    2)
})


test_that("perfectSimulaton does not mine bias change", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 100
  t.vec <- seq(tnum)
  true.change.point <- 40
  true.bias.point <- 60
  before.change.vec <- rnbinom(6*(true.change.point), 30, 0.995)
  p <- 0.05
  after.change.vec <- rnbinom(6*(true.bias.point - true.change.point), 30, 1 - p)
  mean.bias <- 3
  biased.p <- mean.bias*p/(1+(mean.bias-1)*p)
  after.bias.vec <- rnbinom(6*(tnum - true.bias.point), 30, 1 - biased.p)
  count.mat <- Matrix(c(before.change.vec, after.change.vec, after.bias.vec), nrow=6, sparse=T)
  mean.bias.vec <- c(rep(1, true.bias.point),
                     rep(mean.bias, tnum - true.bias.point))
  mean.bias.vec <- mean.bias.vec/mean(mean.bias.vec)
  ## estimation
  nbcf <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                          alpha, beta, r, lambda,
                          p.res, count.res, t.res)
  ## check change points num
  median.change.point.num <- median(unlist(purrr::map(nbcf@sim.change.point.list,
                                                      ~ length(.x))))
  expect_equal(
    median.change.point.num,
    1)
})

test_that("findChangeVariate correctly mine changed variates for each change point", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 700
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  rownames(count.mat) <- as.character(seq(nrow(count.mat)))
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  nbcf <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                          alpha, beta, r, lambda,
                        p.res, count.res, t.res)
  change1.df <- nbcf@change.variate.df %>%
    dplyr::filter(change.point - 400 == min(change.point - 400)) %>%
    dplyr::arrange(dplyr::desc(change.point))
  expect_equal(
    sort(as.numeric(change1.df$var[1:6])),
    sort(seq(6))
  )
    
})


test_that("calculateChangePoints correctly mine changed points for each variate", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 700
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  rownames(count.mat) <- as.character(seq(nrow(count.mat)))
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  change.point.df <- calculateChangePoints(count.mat, t.vec, mean.bias.vec,
                                           alpha, beta, r, lambda,
                                           p.res, count.res, t.res)
  v1.ct.vec <- dplyr::filter(change.point.df, var == "1")$change.point
  expect_lt(
    abs(mean(v1.ct.vec - 400)),
    200
  )
  
  v2.ct.vec <- dplyr::filter(change.point.df, var == "7")$change.point
  expect_lt(
    abs(mean(v2.ct.vec - 700)),
    200
  )
    
})

test_that("calculateChangePoints correctly mine changed points for each variate, using perfectsimulaton", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 700
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  rownames(count.mat) <- as.character(seq(nrow(count.mat)))
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  change.point.df <- calculateChangePoints(count.mat, t.vec, mean.bias.vec,
                                           alpha, beta, r, lambda,
                                           p.res, count.res, t.res, method = "sim")
  v1.ct.vec <- dplyr::filter(change.point.df, var == "1")$change.point
  expect_lt(
    abs(mean(v1.ct.vec - 400)),
    200
  )
  
  v2.ct.vec <- dplyr::filter(change.point.df, var == "7")$change.point
  expect_lt(
    abs(mean(v2.ct.vec - 700)),
    200
  )
    
})

test_that("MAP estimate correctly mine two wchange points", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 700
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  nbcf <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                          alpha, beta, r, lambda,
                          p.res, count.res, t.res)
  ## check distance from true points
  min.dist.from.400 <- min(abs(nbcf@map.change.point - 400))
  expect_lt(
    min.dist.from.400,
    10)
  min.dist.from.700 <- min(abs(nbcf@map.change.point - 700))
  expect_lt(
    min.dist.from.700,
    10)
  ## check change points num
  change.point.num <- length(nbcf@map.change.point)
  expect_equal(
    change.point.num,
    2)
  expect_equal(
    order(nbcf@map.change.point),
    seq(change.point.num)
  )
})



test_that("MAP estimate with fix num correctly mine two wchange points", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 700
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  nbcf <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                          alpha, beta, r, lambda,
                          p.res, count.res, t.res, map.fix.num = 2)
  ## check distance from true points
  min.dist.from.400 <- min(abs(nbcf@map.change.point - 400))
  expect_lt(
    min.dist.from.400,
    10)
  min.dist.from.700 <- min(abs(nbcf@map.change.point - 700))
  expect_lt(
    min.dist.from.700,
    10)
  ## check change points num
  change.point.num <- length(nbcf@map.change.point)
  expect_equal(
    change.point.num,
    2)
  expect_equal(
    order(nbcf@map.change.point),
    seq(change.point.num)
  )
})

test_that("bayes.factor and FDR calculation in calculateChangePoints is work well", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  ## unchanged variates
  c2.vec <- rnbinom(6*(tnum), 30, 0.995)
  c2.count.mat <- Matrix(c2.vec, ncol=tnum, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  rownames(count.mat) <- as.character(seq(nrow(count.mat)))
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  change.point.df <- calculateChangePoints(count.mat, t.vec, mean.bias.vec,
                                           alpha, beta, r, lambda,
                                           p.res, count.res, t.res)
  v1.df <- dplyr::filter(change.point.df, as.numeric(var) <= 6)
  v2.df <- dplyr::filter(change.point.df, as.numeric(var) > 6)
  expect_true(
    min(v1.df$bayes.factor) > max(v2.df$bayes.factor)
  )
  expect_true(
    min(v2.df$FDR) > max(v1.df$FDR)
  )    
})

test_that("bayes.factor of calculateNbcf works well", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 1000
  t.vec <- seq(tnum)
  c1.true.change.point <- 400
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  ## unchanged variates
  c2.vec <- rnbinom(6*(tnum), 30, 0.995)
  c2.count.mat <- Matrix(c2.vec, ncol=tnum, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  rownames(count.mat) <- as.character(seq(nrow(count.mat)))
  mean.bias.vec <- rep(1, length(t.vec))
  ## estimation
  nbcf.obj <- calculateNbcf(count.mat, t.vec, mean.bias.vec,
                                   alpha, beta, r, lambda,
                                   p.res, count.res, t.res)
  expect_true(
    nbcf.obj@bf.list[["1"]] > nbcf.obj@bf.list[["7"]]
  )
})
