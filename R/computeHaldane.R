#' Compute genetic distances (in centiMorgans) by the Haldane method
#'
#'    Compute genetic distances (in centiMorgans) by the Haldane method, and assign a group number to the markers based on which cutpoint group they are in.
#'    NOTE: Recombination frequency values should NOT be divided by 2 before being passed to this function, the function takes care of that.
#' @usage computeHaldane(pathData, cutpoints)
#' @param pathData - TSP solution object, produced by the TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput functions
#' @param cutpoints - cutpoint array, as produced by the cutpoints.R/findcuts function.  The default is NULL, which means that pathData consists of a single chromosome and therefore cutpoints are not needed.
#' @return TSP solution object identical to pathData along with genetic distances and group numbers
#' @export
computeHaldane <- function(pathData, cutpoints=NULL)
{
  # If cutpoints == NULL, create a dummy cutpoint that spans the entirety of the cluster.
  if(is.null(cutpoints))
  {
    cutpoints = matrix(c(1,1,length(pathData),1), 1, 4)
  }

  # Length of arrays.
  arrLen = length(pathData[,1])

  # Create a vector to hold the genetic distances and group numbers.
  haldanecM <- array(0, dim=c(arrLen,2))

  # Label the groups.
  groupNum = 1

  # Start at the first cutpoint.
  haldanecM[cutpoints[1,2],1] = 0
  haldanecM[cutpoints[1,2],2] = groupNum

  # The for loop starts at the marker after the first cutpoint.
  for(i in (cutpoints[1,2]+1):arrLen)
  {
    if(i %in% cutpoints[,2])
    {
      haldanecM[i,1] = 0
      groupNum = groupNum + 1
    }
    else
    {
      nextval = as.numeric(haldanecM[i-1]) + 50 * log( 1 / (1 - as.numeric(pathData[i-1,3])))
      haldanecM[i,1] = nextval
    }
    haldanecM[i,2] = groupNum
  }

  # After we reach the end of the list, start at marker #1 and proceed until the marker just before the first cutpoint.
  # Only do this if the first cutpoint is not marker #1.
  limit = cutpoints[1,2] - 1
  if(limit > 1)
  {
    for(i in 1:limit)
    {
      if(i %in% cutpoints[,1] )
      {
        haldanecM[i,1] = 0
        groupNum = groupNum + 1
      }
      else if(i == 1) # For marker #1, we have to use the last value in the list instead of the i-1 values.
      {
        nextval = as.numeric(haldanecM[arrLen]) + 50 * log( 1 / (1 - as.numeric(pathData[arrLen,3])))
        haldanecM[i,1] = nextval
      }
      else
      {
        nextval = as.numeric(haldanecM[i-1]) + 50 * log( 1 / (1 - as.numeric(pathData[i-1,3])))
        haldanecM[i,1] = nextval
      }
      haldanecM[i,2] = groupNum
    }
  }

  # Add the genetic distance array onto the pathData object.
  colnames(haldanecM) <- c("Haldane (cM)", "group")

  newarray = cbind(pathData, haldanecM)

  return(newarray)
}
