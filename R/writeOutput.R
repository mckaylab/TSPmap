#' Write the TSP solution object to a comma-separated value file.
#'
#'    Write the TSP solution object to a comma-separated value file.
#' @usage writeOutput(filename, TSPsol)
#' @param filename - path and output filename as a string
#' @param TSPsol - TSP solution object, as produced by TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput.
#' @return
#'    None
#' @export
writeOutput <- function(filename, TSPsol)
{
  write.table(TSPsol, file = filename, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
}
