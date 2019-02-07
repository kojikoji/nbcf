##' 
##' Simulated time course negative binomial data
##' 
##' This data contains 6 variate time course negative binomial data.
##' elemnts of t.vec are 0.01, 0.02, ..., 1.0, and includes two change points,
##' where each variate can change the parameter p.
##' 
##' @docType data
##' 
##' @usage data(sim.data)
##' 
##' @example
##' sim.data <- data(sim.data)
##' nbcf.obj <- loadData(sim.data$count.mat, sim.data$t)
##' nbcf.obj <- setupNegativeBinomDict(nbcf.obj)
##' nbcf.obj <- calculateQt(nbcf.obj, lambda= 0.01, alpha=1.0, beta=1.0)
##' nbcf.obj <- perfectSimulation(nbcf.obj)
##' nbcf.obj <- estimateMap(nbcf.obj)
##' nbcf.obj <- detectVariate(nbcf.obj)
##' hist(nbcf.obj@change.point)
##' nbcf.obj@change.variate
##' 
##' @creation
##' library(tidyverse)
##' library(Matrix)
##' change.point <- c(0.2, 0.7)
##' change.var <- list(C1=c(1, 2, 3), C2=c(3, 4, 5, 6))
##' p.list <- list(
##'   I1=c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5),
##'   I2=c(0.3, 0.1, 0.9, 0.5, 0.5, 0.5),
##'   I3=c(0.3, 0.1, 0.2, 0.6, 0.1, 0.8))
##' r.val <- 30
##' generate.sample <- function(p.vec) purrr::map(
##'                                             seq(length(p.vec)),
##'                                             ~ list(label=.x, val=rnbinom(1, 30, p.vec[.x])))
##' sim.df <- tibble(t=seq(100)/100) %>% 
##'   dplyr::mutate(
##'            I=ifelse(
##'              t < 0.2,
##'              "I1",
##'              ifelse(
##'                t < 0.7,
##'                "I2",
##'                "I3"))) %>%
##'   dplyr::mutate(
##'            sample = purrr::map(
##'                              I,
##'                              ~ generate.sample(p.list[[.x]]))) %>%
##'   unnest() %>%
##'   dplyr::mutate(
##'            label = purrr::map(
##'                             sample,
##'                             ~ .x[[1]]),
##'            val = purrr::map(
##'                             sample,
##'                             ~ .x[[2]])) %>%
##'   select(-sample)  %>%
##'   unnest() %>%
##'   tidyr::spread(label, val)
##' sim.data <- list(
##'   t.vec = sim.df$t,
##'   count.mat = Matrix(as.matrix(sim.df[,3:8]), sparse = T))
"sim.data"

