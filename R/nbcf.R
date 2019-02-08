## Nbcf class
## Yasuhiro Kojima

##' Data stroge class which contain the information required for finding change points.
##' 
##' This class is a data storage object for calculating change points and changing variate.
##' You should firstly apply calculateQt and after either perfectSimulation or estimateMap.
##' 
##' @slot count.mat dgCMatrix, whose rows represent time, and columns represent variate.
##' @slot t.vec Numeric vector, whose elements represent observed time corresponding to each row of count matrix.
##' @slot mean.bias.vec Numeric vector, whose elements represent expted bias for mean value of each row.
##' @slot params List, named list of parameters, whose mebers must be "lambda", "r", "alpha" and "beta".
##' @slot lnb.dict LnbDict, Storage for precomputated log negative binomial probability. This is a instance of  \code{\link{LnbDict}}.
##' @slot Numeric vector, Qt probability of Yt:n given there is change point in t-1
##' @slot Numeric vector, sim.change.point change points simulated from perfectSimulation
##' @slot Numeric.vector, map.change.point change points simulated from estimateMap
##' @slot Character vector, change.variate named list of changed variates for each map.change.points
##'
##' @import Matrix

setClass(
  "Nbcf",
  representation(
    count.mat="dgCMatrix",
    t.vec="vector",
    mean.bias.vec="vector",
    params="list",
    lnb.dict.mat="list",
    Qt="vector",
    Pst="matrix",
    sim.change.point = "vector",
    map.change.point = "vector",
    change.variate = "list"
  )
)
