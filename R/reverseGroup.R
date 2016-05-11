#' Reverse the order of markers in a group within a TSP solution object.
#'
#'    Reverse the order of markers in a group within a TSP solution object.
#' @usage reverseGroup(markerData, cutpoints, groupIndex, rfmatrix)
#' @param markerData - TSP solution object, produced by the TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput function
#' @param cutpoints - cutpoint array, produced by cutpoints.R/findcuts function
#' @param groupIndex - index of the group
#' @param rfmatrix - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @return
#'    TSP solution object that is identical to markerData, but with the markers of the selected group in reverse order.
#'    NOTE: this function updates the rf values and cumulative rf values, but does not update genetic distance values (if present).
#' @export
reverseGroup <- function(markerData, cutpoints=NULL, groupIndex=1, rfmatrix=NULL)
{
  # If cutpoints == NULL, create a dummy cutpoint that spans the entirety of the cluster.
  if(is.null(cutpoints))
  {
    cutpoints = matrix(c(1,1,length(markerData[,1]),length(markerData[,1])), 1, 4)
  }

  startpoint = as.numeric(cutpoints[groupIndex,2])
  endpoint = as.numeric(cutpoints[groupIndex,3])
  groupLen = as.numeric(cutpoints[groupIndex,4])
  totalLen = length(markerData[,1])

  for(i in 1:(groupLen/2))
  {
    # Swap the markers at startpoint and endpoint.
    temp = markerData[startpoint,]
    markerData[startpoint,] = markerData[endpoint,]
    markerData[endpoint,] = temp

    # Update start and end points.
    startpoint = startpoint + 1
    if(startpoint > totalLen)
      startpoint = 1
    endpoint = endpoint - 1
    if(endpoint < 1)
      endpoint = totalLen
  }

  # Set firstpoint and secondpoint.
  firstpoint = as.numeric(cutpoints[groupIndex,2])
  if(firstpoint == totalLen)
    secondpoint = 1
  else
    secondpoint = firstpoint + 1

  # Shift rf values.
  for(i in 1:(groupLen-1))
  {
    markerData[firstpoint,3] = as.numeric(markerData[secondpoint,3])

    # Update first and second points.
    firstpoint = firstpoint + 1
    if(firstpoint > totalLen)
      firstpoint = 1
    secondpoint = secondpoint + 1
    if(secondpoint > totalLen)
      secondpoint = 1
  }

  # Look up the new value for the rf value between the marker before startpoint and the marker at startpoint.
  # Set startpoint and prevpoint.
  startpoint = as.numeric(cutpoints[groupIndex,2])
  if(startpoint == 1)
    prevpoint = totalLen
  else
    prevpoint = startpoint - 1

  if(markerData[prevpoint,1] < markerData[startpoint,1])
  {
    lowIndex = as.numeric(markerData[prevpoint,1])
    highIndex = as.numeric(markerData[startpoint,1])
  }
  else
  {
    lowIndex = as.numeric(markerData[startpoint,1])
    highIndex = as.numeric(markerData[prevpoint,1])
  }
  markerData[prevpoint,3] = as.numeric(rfmatrix[highIndex, lowIndex])

  # Also look up the new value for the rf value between the marker at endpoint and the marker after endpoint.
  #Set endpoint and nextpoint.
  endpoint = as.numeric(cutpoints[groupIndex,3])
  if(endpoint == totalLen)
    nextpoint = 1
  else
    nextpoint = endpoint + 1

  if(markerData[endpoint,1] < markerData[nextpoint,1])
  {
    lowIndex = as.numeric(markerData[endpoint,1])
    highIndex = as.numeric(markerData[nextpoint,1])
  }
  else
  {
    lowIndex = as.numeric(markerData[nextpoint,1])
    highIndex = as.numeric(markerData[endpoint,1])
  }
  markerData[endpoint,3] = as.numeric(rfmatrix[highIndex, lowIndex])

  # Update cumulative rf values (column 4).
  # Set firstpoint and secondpoint.
  firstpoint = as.numeric(cutpoints[groupIndex,2])
  if(firstpoint == totalLen)
    secondpoint = 1
  else
    secondpoint = firstpoint + 1

  for(i in 1:groupLen)
  {
    markerData[startpoint,4] = as.numeric(markerData[secondpoint,4]) + as.numeric(markerData[firstpoint,3])

    # Update first and second points.
    firstpoint = firstpoint + 1
    if(firstpoint > totalLen)
      firstpoint = 1
    secondpoint = secondpoint + 1
    if(secondpoint > totalLen)
      secondpoint = 1
  }

  return(markerData)
}

