context("test-calculateChangePoint")

test_that("test calculateCtTransitionProb",{
          lqt <- runif(10, -10, -1)
          lpst <- matrix(runif(10*10, -10, -1), nrow=10)
          lambda <- 0.01
          transition.prob.list <- calculateCtTransitionProb(lqt, lpst, lambda)
          expect_equal(
            length(transition.prob.list),
            10)
          expect_equal(
            sum(transition.prob.list[["2"]]),
            1)
})

test_that("test perfectSimulation",{
          lqt <- runif(10, -10, -1)
          lpst <- matrix(runif(10*10, -10, -1), nrow=10)
          lambda <- 0.5
          change.points.list <- perfectSimulation(lqt, lpst, lambda, sim.iter=5)
          expect_equal(
            length(change.points.list),
            5)
})
