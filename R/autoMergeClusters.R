#' Merge the clusters from autoClusterMST into full linkage groups.
#'
#'  Starting with the clusters produced by autoClusterMST, this function merges them until we are left with full linkage groups.
#' @usage autoMergeClusters(clusterlist, rfmat, LKHexec,
#' numChromosomes, rfCap = 0, minMergeWeight = 0.40, rawdata)
#' @param clusterlist - list of clusters produced by autoClusterMST
#' @param rfmat - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param LKHexec - path and filename of LKH executable - must not contain the tilde character!
#' @param numChromosomes - number of chromosomes
#' @param rfCap - value at which the rf matrix was capped.
#' @param minMergeWeight - this is the minimum weight that will cause two clusters to be merged.  If the min weight between the two groups exceeds this threshold, the two groups will not be merged.
#' @param rawdata original raw marker data that was read into the R session
#' @return Returns the list of linkage groups.
#' @export
autoMergeClusters <- function(clusterlist, rfmat, LKHexec, numChromosomes, rfCap = 0, minMergeWeight = 0.40, rawdata)
{
  # See if the smallest clusters have a clear affinity for another cluster.

  numClusters = length(clusterlist)

  # Define the clusterOrder to hold the cluster indices in order of size.
  clusterSizes = vector()

  # Create a matrix of TSP result values so we don't have to redo any TSP runs that we've already made.  Initialize all values to 100.
  tspvals = matrix(100.0,nrow=numClusters, ncol=numClusters)

  for(i in 1:length(clusterlist))
  {
    clusterSizes = c(clusterSizes, length(clusterlist[[i]]))
  }

  sortedSizes = sort(clusterSizes)

  clusterOrder = order(clusterSizes)

  # Now clusterOrder holds the indices of each cluster in size order.

  # ===========================================================
  # Begin merging groups:

  # Sort the clusters in order of size.  This makes it possible to do the clusters in order with the while loop.
  # If the clusterlist has been produced by the autoGroupMST function they should already be sorted, but we can't be sure that this is the case.

  numClusters = length(clusterlist)

  clusters = list()

  for(i in 1:numClusters)
  {
    clusters = c(clusters, clusterlist[clusterOrder[i]])
  }

  # Initialize testIndex.
  testIndex = 1

  # matchweightlist keeps track of the weights that trigger each merge operation.  This is so that we can inspect these values and see if there were any groups merged which should not have been.
  matchweightlist = vector()

  while(testIndex <= numClusters && numClusters > numChromosomes)
  {
    # This flag indicates whether or not we need to merge the groups at testIndex and smallestIndex at the end of the while loop.
    mergeFlag = FALSE

    # First we find the smallest edge between the cluster at testIndex and all other clusters.  If all clusters have large values for the smallest possible edge weight except for one, we can safely assume that the current cluster should be merged with it.

    # Create a list of the smallest possible edge weights between the current cluster and all others.
    smallestEdges = as.numeric(checkNongroupEdges(rfmat, clusterOrder[testIndex], clusters))

    # The index of the cluster with the smallest weight edge:
    lowestIndex = which(smallestEdges == min(smallestEdges))

    # Sort smallestEdges and check the value of the smallest against the value of the second-smallest.
    sortedSmallestEdges = sort(as.numeric(smallestEdges))

    # Check smallestEdges for the cluster at lowestIndex. This isn't needed and is just a test.
    minweight = min(smallestEdges)

    # If the smallest member of smallestEdges is less than or equal to half the value of the second-smallest member, we have a match.
    if(sortedSmallestEdges[1] <= 0.5*sortedSmallestEdges[2])
    {
      mergeFlag = TRUE
    } else {
      # If that didn't work, we try to pair the cluster at testIndex with all other clusters and run each combination through the TSP solver.  Then we check the largest rf value in the solution, and if one group indicates a very small value we know we can merge the groups.

      # But before doing that, check the smallest edge between the two groups, if it's above the merge threshold then you don't need to actually do the TSP.

      weightlist = vector()

      for(i in 1:numClusters)
      {
        if(i == testIndex)
        {
          testweight = 100
        } else if(tspvals[testIndex,i] == 100){
          # If we already have the result of the TSP run, use it so we don't waste time running the solver.  A value of 100 indicates that we do not yet have the value.
          #cat("\n\n-----------> checking ", i, "\n\n")
          #cat("\n\n-------> testIndex = ", testIndex, "\n")
          cat("number of Clusters = ", numClusters, "\n\n")

          tspvals[i,testIndex] = matchGroups2(clusters[[testIndex]], clusters[[i]], rawdata, rfmat, LKHexec)
          tspvals[testIndex,i] = tspvals[i,testIndex]
        }
      }

      lowestIndex = which(tspvals[testIndex,] == min(tspvals[testIndex,]))

      # It's possible that there is a tie for the minimum weight.
      if(length(lowestIndex > 1))
      {
        lowestIndex = lowestIndex[1]
      }

      # Due to the capping of rf values, we can end up with several weights that are identical, so we need to check that the lowest value is actually smaller than the capped value.
      # Note that this is not the same as verifying that the smallest value is unique: there could be two matching smallest values in the case that a chromosome has been split into 3 clusters.

      minweight = min(tspvals[testIndex,])

      # If minweight is very small, we assume the groups can be merged directly.
      sortedSmallestEdges = sort(as.numeric(tspvals[testIndex,]))

      if(sortedSmallestEdges[1] <= 0.5*sortedSmallestEdges[2])
      {
        mergeFlag = TRUE
      } else if(minweight < minMergeWeight) {
        # Otherwise, we check to see if the cluster at the lowestIndex also shows affinity for the test cluster.
        # Now do the reverse for the group at the smallest rf value in weightlist.
        weightlist = vector()

        for(i in 1:numClusters)
        {
          if(i == lowestIndex)
          {
            testweight = 100
          } else if(tspvals[lowestIndex,i] == 100) {
            tspvals[lowestIndex,i] = matchGroups2(clusters[[lowestIndex]], clusters[[i]],rawdata = rawdata, rfmat, LKHexec)
            tspvals[i,lowestIndex] = tspvals[lowestIndex,i]
          }
        }

        minweight = min(tspvals[lowestIndex,])
        newlowestIndex = which(tspvals[lowestIndex,] == minweight)

        # At this point we check to see if the original cluster (at index testIndex) and the current cluster (at index lowestIndex) each indicate that the other is it's lowest TSP affinity.
        # However, due to the capping of rf values, we can end up with several weights that are identical, so we need to check that the lowest value is actually smaller than the capped value.
        # Note that this is not the same as verifying that the smallest value is unique: there could be two matching smallest values in the case that a chromosome has been split into 3 clusters.

        if(testIndex == newlowestIndex && minweight <= minMergeWeight)
        {
          mergeFlag = TRUE
        }
      }
    }

    # If we found a reason to merge the groups at testIndex and lowestIndex, we now merge them.
    if(mergeFlag)
    {
      print(paste0("merging clusters ", testIndex, " and ", lowestIndex))

      # Merge the smaller cluster (at testIndex) into the larger cluster (at lowestIndex).
      clusters[[lowestIndex]] = c(clusters[[testIndex]], clusters[[lowestIndex]])

      # Delete the smaller cluster
      clusters[[testIndex]] = NULL

      matchweightlist = c(matchweightlist, minweight)

      # Update tspvals matrix.  Reset the values of the row/column of the cluster that was removed, then delete the row/column of the group that was removed (do it in that order).
      tspvals[lowestIndex,] = 100.0
      tspvals[,lowestIndex] = 100.0
      tspvals = tspvals[-testIndex,]
      tspvals = tspvals[,-testIndex]

      testIndex = testIndex - 1
      numClusters = numClusters - 1
    }

    testIndex = testIndex + 1
  }
  return(clusters)
}
