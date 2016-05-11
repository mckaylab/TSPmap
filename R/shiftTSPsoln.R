#' Shift a TSP solution array when the endpoints don't line up with the beginning of the array.
#'
#'    This function returns a list of the smallest edge weight between a given cluster and all other clusters.  This is used to make the decision to merge clusters when the smallest edge weight is below a given threshold.  You will want to re-run the findCuts function after running this method.
#' @usage shiftTSPsoln(TSPsoln, startIndex)
#' @param TSPsoln - array of markers of a TSP solution
#' @param startIndex - StartIndex is the index of the first cutpoint, which we will shift so that it is in position 1.  If startIndex is not specified, the default value of 0 will cause the method to shift the solution array so that it ends at the largest rf value (since this is probably the best cutpoint).
#' @return Array of markers in TSP solution order, but shifted to startIndex.
#' @export
shiftTSPsoln <- function(TSPsoln, startIndex = 0)
{
  # If the startIndex parameter is greater than zero.
  if(startIndex == 0)
  {
    # Find the largest rf value.
    maxrf = max(TSPsoln[,3])
    startIndex = which(TSPsoln[,3] == maxrf)
  }

  # There are situations where startIndex may be a list of more than one item.
  # If that happens, we only want to use the first element.
  if(length(startIndex) > 1)
  {
    startIndex = startIndex[1]
  }

  maxInd = length(TSPsoln[,1])
  for(i in 1:startIndex)
  {
    # Save the first row.
    temp = TSPsoln[1,]

    # Remove the first row.
    TSPsoln = TSPsoln[-1,]

    # Add the first row to the end.
    TSPsoln = rbind(TSPsoln, temp)
  }

  # Clear out row.names before returning.
  row.names(TSPsoln) = NULL

  # Need to recalculate cumulative rf values.
  TSPsoln[1,4] = 0
  for(i in 2:maxInd)
  {
    TSPsoln[i,4] = as.numeric(TSPsoln[(i-1),4]) + as.numeric(TSPsoln[(i-1),3])
  }

  return(TSPsoln)
}
