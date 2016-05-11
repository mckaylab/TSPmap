#' Create minimum spanning tree.
#'
#' This function creates a minimum spanning tree and breaks it into numcuts clusters.
#' @usage createMST(rfmat, numcuts)
#' @param rfmat - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param numcuts - number of clusters you wish to break the MST into.
#' @return List of numcuts+1 marker clusters.
#' @import igraph
#' @export
createMST <- function(rfmat, numcuts)
{
  # before doing anything, we have to disallow the diagonal elements since these will all be zero.
  for(i in 1:length(rfmat[,1]))
  {
    rfmat[i,i] = 100
  }

  g1 <- graph.adjacency(rfmat, weighted = T, mode = "undirected")

  mymst <- minimum.spanning.tree(g1)

  # print weights in descending order.
  print(sort((E(mymst)$weight), decreasing=TRUE)[1:numcuts])

  weightcutoff = sort((E(mymst)$weight), decreasing=TRUE)[numcuts]

  mymstcopy = delete.edges(mymst, which(E(mymst)$weight >= weightcutoff))

  # Get cluster info.
  clusterinfo = clusters(mymstcopy)
  clustersizes = clusterinfo$csize
  sortedclustersizes = sort(clustersizes)

  return(clusterinfo)
}
