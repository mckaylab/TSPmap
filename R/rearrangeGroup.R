#' Merge two clusters.
#'
#'    Merge two clusters.  This is done after finding the results of matchGroups and knowing which two groups are actually part of the same chromosome.=
#' @usage rearrangeGroup(group1Index, group2Index, cutlist, pathData)
#' @param group1Index - index of the first group
#' @param group2Index - index of the second group
#' @param cutlist - cutlist array, produced by cutpoints.R/findcuts function
#' @param pathData - TSP solution object, produced by TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput
#' @return
#'    List of the new cutpointArray and TSP solution, after the two groups have been merged.
#'    Get the return values like this:
#'    result = rearrangeGroup(1, 2, cutpointArray, TSPsolution)
#'    newCutpointArray = result$v1
#'    newTSPsolution = result$v2
#' @export
rearrangeGroup <- function(group1Index, group2Index, cutlist, pathData)
{
  # Find the smaller group.
  if(cutlist[group1Index,4] < cutlist[group2Index,4])
  {
    smallIndex = group1Index
    largeIndex = group2Index
  }
  else
  {
    smallIndex = group2Index
    largeIndex = group1Index
  }

  # Find the group with the lower index.
  if(group1Index < group2Index)
  {
    lowIndex = group1Index
    highIndex = group2Index
  }
  else
  {
    lowIndex = group2Index
    highIndex = group1Index
  }

  # Combine the smaller group with the larger.
  smallSize = cutlist[smallIndex,4]

  # Increase the size of the larger group.
  cutlist[largeIndex,4] = cutlist[largeIndex,4] + smallSize

  # If the smaller group also has the lower index:
  if(lowIndex == smallIndex)
  {
    # Decrease the starting point of the larger group.
    cutlist[largeIndex,2] = cutlist[largeIndex,2] - smallSize

    # Decrease the start/end point of groups below highIndex.
    # Only do this if highIndex is more than 1 greater than smallIndex.
    if(lowIndex < (highIndex - 1))
    {
      for(i in (highIndex-1):(lowIndex+1))
      {
        cutlist[i,2] = cutlist[i,2] - smallSize
        cutlist[i,3] = cutlist[i,3] - smallSize
      }
    }
  }
  else  # -=-=-=-= Otherwise the smaller group has the higher index, so do the inverse. -=-=-=-=
  {
    # Increase the ending point of the larger group.
    cutlist[largeIndex,3] = cutlist[largeIndex,3] + smallSize

    # Increase the start/end points of groups above highIndex, but only if smallIndex is more than 1 greater than highIndex.
    if(highIndex > (lowIndex + 1))
    {
      for(i in (lowIndex+1):(highIndex-1))
      {
        cutlist[i,2] = cutlist[i,2] + smallSize
        cutlist[i,3] = cutlist[i,3] + smallSize
      }
    }
  }

  # Set the starting markerID of largeIndex to be equal to the starting markerID of smallIndex.
  cutlist[largeIndex,1] = cutlist[smallIndex,1]


  # Move all markers from the small group to the beginning of the large group in markerData.
  moveGroup = pathData[cutlist[smallIndex,2]:cutlist[smallIndex,3],]

  newpathData = pathData[c(1:(cutlist[smallIndex,2]-1),(cutlist[smallIndex,3]+1):length(pathData[,1])),]

  pathSegment1 = newpathData[1:(cutlist[largeIndex,2]-1),]
  pathSegment2 = newpathData[(cutlist[largeIndex,2]):length(newpathData[,1]),]

  finalPath = rbind(pathSegment1, moveGroup, pathSegment2)

  # Remove the smaller group from the culist.
  cutlist = cutlist[-smallIndex,]

  return(list(v1=cutlist, v2=finalPath))
}
