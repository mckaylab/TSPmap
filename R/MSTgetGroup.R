#' Get the marker IDs for a particular cluster
#'
#'    This function returns the list of marker IDs for a given cluster.
#' @usage MSTgetGroup(clusterinfo,clusterIndex)
#' @param clusterinfo list of clusters, as produced by createMST function.
#' @param clusterIndex  the index of the cluster.
#' @return Marker IDs of the selected cluster.
#' @export
MSTgetGroup <- function(clusterinfo, clusterIndex)
{
  cluster = which(clusterinfo$membership == clusterIndex)
  return(cluster)
}
