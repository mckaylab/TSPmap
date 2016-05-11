#' Process LKH solution file.
#'
#'    Process LKH solution file.
#' @usage processLKHOutput(filename, markerdata, rfmat)
#' @param filename - path and name of LKH solution file, which is returned from the function callLKH (this file has a .LKH extension)
#' @param markerdata - main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function
#' @param rfmat - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @return
#'    TSP solution object, containing:
#'      1. the index of the markers in the main marker data object
#'      2. marker name
#'      3. rf value between each marker and the next
#'      4. cumulative rf value
#' @export
processLKHOutput <- function(filename, markerdata, rfmat)
{
  # The LKH output file has a header that needs to be removed after we read it in.
  # The path terminates with a "-1".

  # Skip the first 6 lines since these are the header.
  # Only read in length(markerdata) values since this is the length of the path.  This avoids the problem of reading in the terminal -1 and "EOF".
  # Add 1 to the value of n because of the dummy marker.
  path <- scan(file=filename, skip=6, n = length(rfmat[,1])+1)

  # Create a list of marker names.
  markernames = array()

  # We want to remove the dummy marker from path.  This marker always has the largest index value, which is 1 greater than the length of markerdata.
  path = path[-which(path == (length(markerdata)+1))]

  for (i in 1:length(path))
  {
    # This grabs the marker names from markerdata.
    markernames[i] = colnames(markerdata)[path[i]]
  }

  # Get rf values and cumulative rf from addRfs function.
  rflist = addRfs(path, rfmat)

  # Bind them all together.
  markeroutput = cbind(path, markernames, rflist)

  # Set column names.
  colnames(markeroutput) <- c("index", "markername", "rf", "cumulative rf")

  # Shift the solution to the endpoint of the chromosome.
  returnlist = shiftTSPsoln(markeroutput)

  return(returnlist)
}
