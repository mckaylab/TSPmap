#'  This function pares down a given rf matrix by keeping or removing the markers that are not in the provided subtour.
#'
#'    This function pares down a given rf matrix by keeping or removing the markers that are not in the provided subtour.
#'    This is used when we want to pass a subset of the markers into a TSP solver.
#' @usage pareMatrix(rfmatrix, markerIndexList, keep, clusters)
#' @param rfmatrix - matrix of recombination frequency values produced by the rfmatrix.R/computeRFmat function
#' @param markerIndexList - list of the indices of the markers that you want to be included in the pared down matrix.  NOTE: this list is most easily obtained by using the cutpoint array that is produced by the cutpoints.R/findcuts function.
#' @param keep (OPTIONAL) - this flag indicates whether you want to keep the markers that are in the markerIndexList, or discard them.
#' @param clusters (OPTIONAL) - this flag indicates whether you are passing in a cluster of marker ids or a markerlist.
#' @return submatrix of rfmatrix that only includes or excludes the markers given in markerIndexList, depending on the value of the keep flag.
#' @export
pareMatrix <- function(rfmatrix, markerIndexList, keep = TRUE, clusters = FALSE)
{
  if(clusters)
  {
    markers = unlist(markerIndexList)
  } else {
    markers = as.numeric(markerIndexList)
  }

  if(keep)
  {
    returnMat = rfmatrix[markers,markers]
  } else {
    returnMat = rfmatrix[-c(markers),-c(markers)]
  }

  return(returnMat)
}


