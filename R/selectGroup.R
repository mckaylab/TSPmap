#' Return the list of marker indices that comprise the given group.
#'
#'    Return the list of marker indices that comprise the given group.
#'    This is an internal method that is used in the matchGroups function.
#'    GROUP SHOULD BE SHIFTED TO BEGIN AT A CUTPOINT BEFORE CALLING THIS???
#' @usage selectGroup(groupIndex, markerData, cutList)
#' @param groupIndex - index of the group to select
#' @param markerData - main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function
#' @param cutList - cutpoint array, produced by cutpoints.R/findcuts function
#' @return List of indices of the markers in the given group.
#' @export
selectGroup <- function(groupIndex, markerData, cutList)
{
  maxInd = length(cutList[,1])

  # For the last group, we may have to accout for wraparound.
  # If there is wraparound, the end index will be smaller than the beginning index.
  if(cutList[groupIndex,3] < cutList[groupIndex,2])
  {
    firstGroup = as.numeric(markerData[cutList[groupIndex,2]:length(markerData[,1])])
    secondGroup = as.numeric(markerData[1:cutList[groupIndex,3]])

    return(c(firstGroup, secondGroup))
  }
  else
  {
    return(as.numeric(markerData[cutList[groupIndex,2]:cutList[groupIndex,3]]))
  }
}
