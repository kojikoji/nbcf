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
##' @slot lpst Numeric matrix, log probability of Ys:t
##' @slot lqt Numeric vector, log probability of Yt:n given there is change point in t-1
##' @slot lpst.list Numeric matrix, list of \code{lpst} for each variate
##' @slot lqt.list Numeric matrix, list of \code{lqt} for each variate
##' @slot bf.list Numeric matrix, list of bayes factor for each variate. The larger values indicate the variate changed in time course.
##' @slot used.vars Character vector, variates used for change points calculation
##' @slot sim.change.point.list List, each elements contains change points simulated from perfectSimulation
##' @slot map.change.point Numeric.vector, change points simulated from estimateMap
##' @slot change.variate.df tbl, data frame of changed variates. There are columns for change points, variate and etropy of change points. Entropy is expected to be lower for moredistinct change points.
##'
##' @import Matrix
##' @import tibble

setClass(
  "Nbcf",
  representation(
    count.mat="dgCMatrix",
    t.vec="vector",
    mean.bias.vec="vector",
    params="list",
    lpst="matrix",
    lqt="vector",
    lpst.list="list",
    lqt.list="list",
    bf.list="list",
    used.vars="vector",
    sim.change.point.list = "list",
    map.change.point = "vector",
    change.variate.df = "tbl"
  )
)
