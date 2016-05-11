#' Determine cutpoints in TSP solution object.
#'
#'    Determine cutpoints in TSP solution object.
#' @usage findcuts(tourData, numCuts, alt)
#' @param tourData - TSP solution object, produced by the TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput function
#' @param numCuts - number of cutpoints to find
#' @param alt - alternate mode, where we instead find the cutpoints that are above a given threshold (specified by the numCuts parameter).
#' @return
#'    array of cutpoint info, containing:
#'      1. marker index of the marker at the cutpoint
#'      2. starting index of the group
#'      3. ending index of the group
#'      4. size of the group
#' @export
findcuts <- function(tourData, numCuts, alt = FALSE)
{
  # Length of TSP tour.
  maxindex = length(tourData[,1])

  if(alt)
  {
    # Create a list of the marker indices that have weights above the threshold specified in numCuts.
    topMarkers = tourData[which(as.numeric(tourData[,3]) >= numCuts),2]
    # Now the number of cuts is the length of topMarkers.  This is used later in this script.
    numCuts = length(topMarkers)
  } else {
    # Create a list of the marker indices that correspond to the largest weights in the tourData array.
    topMarkers = tourData[order(-as.numeric(tourData[,3])),2][1:numCuts]
  }

  cutList = array()

  for(i in 1:length(topMarkers))
  {
    # Add 1 here since the cut point is the second marker of the two that have a large weight between them.
    cutIndex = which(tourData[,2] == topMarkers[i])+1

    # If there are no valid cutpoint, we set cutIndex to 1.
    if(length(cutIndex) == 0)
    {
      cutIndex = 1
    }

    if(cutIndex > length(tourData[,1]))
      cutIndex = 1
    cutList[i] = as.integer(tourData[cutIndex,1])

    # Index of the marker in the sorted RF matrix.
    cutloc = as.integer(which(tourData[,2] == topMarkers[i]))
  }

  # Find the indices of the cutpoints in tourData.
  cutind = array()
  for(i in 1:numCuts)
  {
    cutind[i] = which(tourData[,1] == cutList[i])
  }

  # Sort the cutpoint indices.
  cutind = sort(cutind)

  # Sort the marker cutpoints in the order that they appear in tourData.
  sortedInd = array()
  for(i in 1:numCuts)
  {
    sortedInd[i] = tourData[cutind[i],1]
  }

  # Find the sizes of the groups.
  groupSizes = array()

  # If we only have 1 cutpoint.
  if(numCuts == 1)
  {
    groupSizes[i] = maxindex

    # Ending indices of the groups.
    EndInd = (cutind[1]+maxindex-1) %% maxindex
    if(EndInd == 0)
    {
      EndInd = maxindex
    }
  } else {
    # If we have multiple cutpoints.
    for(i in 2:numCuts)
    {
      groupSizes[i-1] = cutind[i] - cutind[i-1]
    }

    # Set the last value of group sizes by calculating the size from the last cutpoint to the first.
    groupSizes[i] = cutind[1] + maxindex - cutind[i]

    # Ending indices of the groups.
    EndInd = cutind[c(2:length(cutind),1)]-1
    # For ending indices which have the value of 0, set them to instead have the value of maxInd.
    EndInd[which(EndInd == 0)] = maxindex
  }

  returnList = cbind(as.numeric(sortedInd), as.numeric(cutind), as.numeric(EndInd), as.numeric(groupSizes))

  colnames(returnList) <- c("first marker ID", "start index", "end index", "size")

  return(returnList)
}
