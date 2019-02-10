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
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec, r, alpha, beta, count.res, p.res)
  ## t.grids will seprate observation into t.res bins
  t.grids <- as.integer(seq(0, length(t.vec), length.out = t.res+1))
  ## count.mat is ordered based on t.vec 
  count.mat <- count.mat[, order(t.vec)]
  lpst <- calculateLpst(lnb.dict, count.mat, t.grids)
  lqt <- calculateLqt(lpst, lambda)
  sim.change.point <- perfectSimulation(lqt, lpst, lambda)
})


test_that("perfectSimulaton correctly mine two wchange points", {
  alpha=2.0
  beta=2.0
  r=30
  lambda=1.0e-2
  p.res=1000
  count.res=500
  t.res=100
  tnum <- 100
  t.vec <- seq(tnum)
  c1.true.change.point <- 40
  c1.before.change.vec <- rnbinom(6*(c1.true.change.point), 30, 0.995)
  c1.after.change.vec <- rnbinom(6*(tnum - c1.true.change.point), 30, 0.95)
  c1.count.mat <- Matrix(c(c1.before.change.vec, c1.after.change.vec), nrow=6, sparse=T)
  c2.true.change.point <- 70
  c2.before.change.vec <- rnbinom(6*(c2.true.change.point), 30, 0.995)
  c2.after.change.vec <- rnbinom(6*(tnum - c2.true.change.point), 30, 0.95)
  c2.count.mat <- Matrix(c(c2.before.change.vec, c2.after.change.vec), nrow=6, sparse=T)
  count.mat <- rbind(c1.count.mat, c2.count.mat)
  mean.bias.vec <- rep(1, length(t.vec))
  lnb.dict <- setupLnbDict(max(count.mat), mean.bias.vec, r, alpha, beta, count.res, p.res)
  ## t.grids will seprate observation into t.res bins
  t.grids <- as.integer(seq(0, length(t.vec), length.out = t.res+1))
  ## count.mat is ordered based on t.vec 
  count.mat <- count.mat[, order(t.vec)]
  lpst <- calculateLpst(lnb.dict, count.mat, t.grids)
  lqt <- calculateLqt(lpst, lambda)
  sim.change.point <- perfectSimulation(lqt, lpst, lambda)
})
