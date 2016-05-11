#' Remove redundant markers.
#'
#'    This function removes redundant markers from a TSP solution list.
#' @usage removeMarkers(TSPsolution, redundantMarkerList, rfmat)
#' @param TSPsolution - TSP solution object, produced by the TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput functions
#' @param redundantMarkerList - list of redundant markers
#' @param rfmat - rf matrix
#' @return TSP solution with redundant markers removed.
#' @export
removeMarkers <- function(TSPsolution, redundantMarkerList, rfmat)
{
  TSPsolution = TSPsolution[-redundantMarkerList,]

  print(length(TSPsolution[,1]))

  # Get new rf and cumulative rf values.
  rfvals = addRfs(TSPsolution[,1], rfmat)

  returnlist = cbind(TSPsolution[,(1:2)],rfvals)

  return(returnlist)
}
