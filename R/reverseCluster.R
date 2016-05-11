#' Reverse the order of markers in a group within a cluster.
#'
#'    Reverse the order of markers in a group within a  cluster
#' @usage reverseCluster(cluster, rfmat)
#' @param cluster cliuster object
#' @param rfmat recombination frequency matrix, produced by the computeRFmat function
#' @return
#'    cluster that is identical to cluster object, but with the markers of the selected  in reverse order.
#'    NOTE: this function updates the rf values and cumulative rf values, but does not update genetic distance values (if present).
#' @export
reverseCluster <- function(cluster, rfmat)
{
  clusterLen = length(cluster[,1])

  for(i in 1:(clusterLen/2))
  {
    # Swap the markers at i and (clusterLen - i + 1)
    temp = cluster[i,]
    cluster[i,] = cluster[(clusterLen - i + 1),]
    cluster[(clusterLen - i + 1),] = temp
  }

  # Update the rf values in columns 3 and 4.

  cumRf = 0

  cluster[1,4] = cumRf

  for(i in 1:(clusterLen-1))
  {
    cluster[i,3] = rfmat[as.integer(cluster[i,1]), as.integer(cluster[(i+1),1])]
    cumRf = cumRf + as.numeric(cluster[i,3])
    cluster[(i+1),4] = cumRf
  }

  # For the last element:
  cluster[clusterLen,3] = rfmat[as.integer(cluster[clusterLen,1]), as.integer(cluster[1,1])]
  cumRf = cumRf + as.integer(cluster[clusterLen,3])

  return(cluster)
}
