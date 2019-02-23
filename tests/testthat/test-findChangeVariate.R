context("test-findChangeVariate")

test_that("calculateEntropyTwoCt works without error",{
          lqt <- runif(10, -10, -1)
          lpst <- matrix(runif(10*10, -10, -1), nrow=10)
          lambda <- 0.5
          entropy <- calculateEntropyTwoCt(lqt, lpst, lambda)
          expect_gte(
              entropy,
              0
          )
})

test_that("findChangeVariate returns right shape data frame",{
    lpst <- matrix(runif(10*10, -10, -1), nrow=10)
    lpst.list <- list(v1 = lpst, v2 = lpst, v3 = lpst)
    ct.vec <- c(3, 8)
    used.vars <- c("v1", "v2")
    change.var.df <- findChangeVariate(lpst.list, used.vars, ct.vec)
    expect_equal(
        nrow(change.var.df),
        2 * 2)
    expect_equal(
        ncol(change.var.df),
        3)

})
