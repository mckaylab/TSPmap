#' This is an alternate version of the matchGroups function above.  This version works with groups from the GroupList.
#'
#'    Combine two cutpoint groups and pass them into LKH and return the value of the second largest rf value in the result.
#'    This is used to pair small subgroups with a larger group by combining the small group with all larger groups and seeing which results in the smallest return value.
#'    The reason we take the second largest rf value in the TSP result is because the largest rf value represents the completion of the TSP circuit, so even for a single chromosome this will be a large value.
#' @usage matchGroups2(group1, group2 , rawdata, rfmatrix, execPath)
#' @param group1 index of the first cluster
#' @param group2 - index of the second cluster
#' @param rawdata - main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function
#' @param rfmatrix - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param execPath - path and name of LKH executable
#' @return Maximum rf value in the TSP path of the two clusters combined.
#' @export
matchGroups2 <- function(group1, group2, rawdata, rfmatrix, execPath)
{
  subgroups = sort(as.numeric(c(group1, group2)))

  # Find the submatrix of the rf matrix that includes all these markers.
  submatrix = pareMatrix(rfmatrix, subgroups)

  # Also select the corresponding subgroup of the raw marker data object.
  subrawdata = rawdata[subgroups]

  # NOTE: tempdir() needs to have an extra "/" added to it.
  tempdirectory = paste0(getwd(),"/")

  # Create the TSP matrix file.
  subtspmat = createTSPFile(subrawdata,submatrix,"temp", tempdirectory)

  # Call LKH.
  callLKH(execPath, paste0(tempdirectory,'./temp.tsp'), tempdirectory)

  # Process the LKH result.
  SubMarkerlist = processLKHOutput(paste0(tempdirectory,"./temp.LKH"), subrawdata, submatrix)

  # To find the maxweight, we want to ignore the largest value because this is where the end of the chromosome wraps around to the beginning, so it's not a valid value.
  weights = sort(SubMarkerlist[,3], decreasing = TRUE)
  maxweight = weights[2]

  return(maxweight)
}


