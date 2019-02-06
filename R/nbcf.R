## nbcf class
## Yasuhiro Kojima

##' Data stroge class which contain the information required for finding change points.
##' 
##' This class is a data storage object for calculating change points and changing variate.
##' You should firstly apply calculateQt and after either perfectSimulation or estimateMap.
##' 
##' @slot count.mat Object of class \code{"ANY"} This is a countmatrix, whose rows represent time, and columns represent variate.
##' @slot t.vec This is a numeric vector, whose elements represent observed time corresponding to each row of count matrix.
##' @slot params This is a named list of parameters, whose mebers must be "lambda", "r", "alpha" and "beta".
##' @slot Qt probability of Yt:n given there is change point in t-1
##' @slot nb.dict.mat Storage for precomputated negative binomial probability.
##' @slot sim.change.point change points simulated from perfectSimulation
##' @slot map.change.point change points simulated from estimateMap
##' @slot change.variate named list of changed variates for each map.change.points


nbcf <- setClass(
  "nbcf",
  slots = c(
    count.mat="ANY",
    t.vec="vector",
    params="list",
    Qt="list",
    nb.dict.mat="ANY",
    sim.change.point = "vector",
    map.change.point = "vector",
    change.variate = "list"
  )
)
