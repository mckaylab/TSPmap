#' Find the smallest edge weight between the target cluster and all other clusters.
#'
#'    This function returns a list of the smallest edge weight between a given cluster and all other clusters.  This is used to make the decision to merge clusters when the smallest edge weight is below a given threshold.
#' @usage checkNongroupEdges(rfmatrix, targetIndex, groupList)
#' @param rfmatrix - rf matrix
#' @param targetIndex - index of the target cluster
#' @param groupList - list of all clusters
#' @return List of minimum edge weights between the cluster at targetIndex and all other clusters.
#' @export
checkNongroupEdges = function(rfmatrix, targetIndex, groupList)
{
  targetGroup = as.numeric(groupList[[targetIndex]])

  minList = list()
  for(i in 1:length(groupList))
  {
    if(i == targetIndex)
    {
      minList = c(minList, as.numeric(100))
    } else {
      testGroup = as.numeric(groupList[[i]])

      minList = c(minList, as.numeric(min(rfmatrix[c(targetGroup), c(testGroup)])))
    }
  }

  return(minList)
}
