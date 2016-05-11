#' Identify redundant markers in the solution object that have rf = 0 so that they can be removed.
#'
#'    This function creates a list of markers in the solution marker list which have 0 recombination frequency with their neighbors, so that they can be removed.
#'    NOTE: after running this function, the redundant marker list need to be removed from the solution marker list, the rf matrix, and the raw data object.  This can be done with the removeMarkers function, the pareMatrix function, and the removedups function, respectively.
#' @usage findZeroes(TSPsolution, rfmat)
#' @param TSPsolution - TSP solution object, produced by the TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput functions
#' @param rfmat recombination frequency matrix
#' @return List of marker indices that can be removed.
#' @export
findZeroes <- function(TSPsolution, rfmat)
{
  # Find the indices of all markers in TSPsolution that have rf = 0 between them and their successor.
  zeroindices = which(TSPsolution[,3] == 0)

  # For each of these, if the next list item also has rf = 0, the next marker can be removed because it adds nothing to the solution.
  # If the next list item has rf > 0, choose which one to keep by comparing the previous rf value to the next.  If the previous rf value is smaller, keep the current marker and discard the next.

  redundantlist = array()

  for(i in 1:length(zeroindices))
  {
    ind = zeroindices[i]

    if(zeroindices[i] != -1)
    {
      # j keeps track of which marker is at the end of the redundant block.
      j = 1
      while(TSPsolution[(ind+j), 3] == 0)
      {
        # Add ind+j to redundanlist.
        redundantlist[length(redundantlist) + 1] = (ind+j)

        # Also mark zeroindices[i+j] to -1 so it will be ignored.
        zeroindices[i+j] = -1

        # Increment j.
        j = j + 1
      }

      # Now that all the interior redundant markers have been removed, choose to keep the marker that has less missing calls. To do this, look at the rf values between the two similar markers and the markers that they border.  The lower total value is the one to keep.

      # Find the total rf value between the marker at ind and the markers before and after.
      tot1 = rfmat[as.integer(TSPsolution[ind,1]),as.integer(TSPsolution[ind-1,1])]+ rfmat[as.integer(TSPsolution[ind,1]),as.integer(TSPsolution[ind+j+1,1])]

      # Find the total rf value between the marker at ind+j and the markers before and after.
      tot2 = rfmat[as.integer(TSPsolution[ind+j,1]),as.integer(TSPsolution[ind-1,1])] + rfmat[as.integer(TSPsolution[ind+j,1]),as.integer(TSPsolution[ind+j+1,1])]

      if(tot1 > tot2)
      {
        redundantlist[length(redundantlist) + 1] = ind
      } else {
        redundantlist[length(redundantlist) + 1] = ind+j
      }
    }
  }

  # Remove the first element of redundantlist, since it is "NA".
  return(redundantlist[-1])
}
