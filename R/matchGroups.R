#' Combine two cutpoint groups and pass them into LKH and return the value of the second largest rf value in the result.
#'
#'    Combine two cutpoint groups and pass them into LKH and return the value of the second largest rf value in the result.
#'    This is used to pair small subgroups with a larger group by combining the small group with all larger groups and seeing which results in the smallest return value.
#'    The reason we take the second largest rf value in the TSP result is because the largest rf value represents the completion of the TSP circuit, so even for a single chromosome this will be a large value.
#' @usage matchGroups(ind_group1,ind_group2, markerData, cutList, rfmatrix, execPath)
#' @param ind_group1 - index of the first group
#' @param ind_group2 - index of the second group
#' @param markerData - main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function
#' @param cutList - cutpoint array, produced by cutpoints.R/findcuts function
#' @param rfmatrix - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param execPath - path to LKH executable
#' @return Maximum rf value in the TSP path of the two clusters combined.
#' @export
matchGroups <- function(ind_group1, ind_group2, markerData, cutList, rfmatrix, execPath)
{
  # Select the markers that make up the groups.
  groupA = selectGroup(ind_group1, markerData, cutList)
  groupB = selectGroup(ind_group2, markerData, cutList)

  subgroups = sort(c(groupA, groupB))

  # Find the submatrix of the rf matrix that includes all these markers.
  submatrix = pareMatrix(rfmatrix, subgroups)

  # Also select the corresponding subgroup of the raw marker data object.
  subrawdata = markerData[c(groupA, groupB)]

  # Create the TSP matrix file.
  subtspmat = createTSPFile(subrawdata,submatrix,"temp", "./")

  # Call LKH.
  callLKH(execPath, './temp.tsp', './')

  # Process the LKH result.
  SubMarkerlist = processLKHOutput("./temp.LKH", subrawdata, submatrix)

  # To find the maxweight, we want to ignore the largest value because this is where the end of the chromosome wraps around to the beginning, so it's not a valid value.
  maxweight = max(SubMarkerlist[,3])

  # Remove max value from Submarkerlist.
  maxloc = which(SubMarkerlist[,3] == maxweight)
  SubMarkerlist = SubMarkerlist[-maxloc,]

  # Find new max weight.
  maxweight = max(SubMarkerlist[,3])

  # Remove the temp files.
  if(file.exists("./temp.LKH"))
    file.remove("./temp.LKH")
  if(file.exists("./temp.tsp"))
    file.remove("./temp.tsp")
  if(file.exists("./temp.par"))
    file.remove("./temp.par")

  return(maxweight)
}
