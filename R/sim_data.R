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
##' @examples
##' sim.data <- data(sim.data)
##' nbcf.obj <- loadData(sim.data$count.mat, sim.data$t)
##' nbcf.obj <- setupNegativeBinomDict(nbcf.obj)
##' nbcf.obj <- calculateQt(nbcf.obj, lambda= 0.01, alpha=1.0, beta=1.0)
##' nbcf.obj <- perfectSimulation(nbcf.obj)
##' nbcf.obj <- estimateMap(nbcf.obj)
##' nbcf.obj <- detectVariate(nbcf.obj)
##' hist(nbcf.obj@change.point)
##' nbcf.obj@change.variate
"sim.data"
