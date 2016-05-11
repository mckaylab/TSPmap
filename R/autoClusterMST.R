#' Break the markers into clusters.
#'
#'  split the MST into a number of clusters equal to 1.5 times the number of chromosomes.  A cluster, in this case, must be at least numMarkers/(2*numChromosomes) in size, so that we don't count exceedingly small groups toward the total number.
#' @usage autoClusterMST(rawdata, rfmat, numChromosomes,
#' LKHexec, internalRfThreshold)
#' @param rawdata - raw marker data
#' @param rfmat - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param numChromosomes - number of chromosomes
#' @param LKHexec - path and filename of LKH executable - must not contain the tilde character!
#' @param internalRfThreshold - if a cluster contains an rf value this large, it probably contains markers from two different chromosomes.  The clusters will be split until no rf values exceeding this threshold exist within any cluster.
#' @return Returns the list of clusters of marker IDs.
#' @export
autoClusterMST <- function(rawdata, rfmat, numChromosomes, LKHexec, internalRfThreshold = 0.30)
{
  # Number of markers.
  numMarkers = length(rfmat[1,])

  # Threshold for defining a cluster is numMarkers/(2*numChromosomes)
  clusterThreshold = as.integer(numMarkers/(2*numChromosomes))

  # Start by breaking into 1/5*numChromosome cluster.
  rawNumClusters = as.integer(1.5*numChromosomes)

  numClusters = 0

  # The while loop will run until we either have enough clusters that meet the size criterion, or we have 2*numChromosomes clusters.  This prevents the function from running endlessly.
  while(numClusters < numChromosomes && rawNumClusters < 2*numChromosomes)
  {
    MSTinfo = createMST(rfmat, rawNumClusters-1)

    # result$csize is the vector of cluster sizes.
    numClusters = length(which(MSTinfo$csize >= clusterThreshold))

    # Keep increasing the number of groups until we have 1.5*numChromosomes.
    rawNumClusters = rawNumClusters + 1
  }

  # Now run each group through LKH and check internalRfThreshold.
  # If any cluster has an internal RF value larger than the threshold, break it apart.  Each of the two new clusters needs to be run through the internal RF test again.

  numClusters = length(MSTinfo$csize)

  # Put the clusters into a list.
  i = 1
  clusterlist = list(MSTgetGroup(MSTinfo, i))
  for(i in 2:numClusters)
  {
    clusterlist = c(clusterlist, list(MSTgetGroup(MSTinfo, i)))
  }

  currIndex = 1

  while(currIndex <= numClusters)
  {
    # It's posible that there may be a group with only a single marker in it.
    if(length(clusterlist[[currIndex]]) > 1)
    {
      # Run each clusterlist[i] through LKH and put the results in a list.
      tspsoln = MSTcallTSP(clusterlist[currIndex], rawdata, rfmat, LKHexec)

      # Print out the second largest weight in the tour, as an indicator of whether or not the group contains two chromosomes.
      topval = as.numeric(sort(tspsoln[,3], decreasing = TRUE)[2])

      if(topval > internalRfThreshold)
      {
        tempCluster = shiftTSPsoln(tspsoln)

        cutPosition = which(as.numeric(tempCluster[,3]) == topval)

        newCluster1 = list(as.numeric(tempCluster[1:cutPosition]))
        newCluster2 = list(as.numeric(tempCluster[(cutPosition+1):length(tempCluster[,1])]))

        # Delete the old cluster from clusterlist.
        clusterlist[[currIndex]] = NULL
        # Add the two new clusters.
        clusterlist = c(clusterlist, newCluster1)
        clusterlist = c(clusterlist, newCluster2)

        currIndex = currIndex - 1
        numClusters = numClusters + 1
      }
    }

    # Remove tspsoln at the end of each iteration to prevent problems with array length.
    if(exists("tspsoln"))
    {
      rm(tspsoln)
    }
    # Increase the value of currIndex.
    currIndex = currIndex + 1
  }

  # Sort clusters in order of size before returning.

  # Define the clusterOrder to hold the cluster indices in order of size.
  clusterSizes = vector()

  for(i in 1:length(clusterlist))
  {
    clusterSizes = c(clusterSizes, length(clusterlist[[i]]))
  }

  sortedSizes = sort(clusterSizes)

  clusterOrder = vector()
  for(i in 1:length(sortedSizes))
  {
    clusterOrder[i] = which(clusterSizes == sortedSizes[i])
    # If there are two clusters with the same size, we will get a duplicate here.  To avoid this, set the clusterSizes element to zero after we've used it.
    clusterSizes[clusterOrder[i]] = 0
  }

  # Rearrange clusterlist so it goes in order of increasing size.
  # This makes it possible to do the clusters in order with the while loop.
  numClusters = length(clusterlist)

  clusters = list()

  for(i in 1:numClusters)
  {
    clusters = c(clusters, clusterlist[clusterOrder[i]])
  }

  # Return the ordered clusters.
  return(clusters)
}
