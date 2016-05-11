#'Process Concorde solution file
#'
#'    Process Concorde solution file.
#' @usage processConcordeOutput(filename, markerdata, rfmat)
#' @param filename path and name of Concode solution file, which is returned from the function callConcorde (this file has a .sol extension)
#' @param markerdata main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function
#' @param rfmat recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @return
#'    TSP solution object, containing:
#'      1. the index of the markers in the main marker data object
#'      2. marker name
#'      3. rf value between each marker and the next
#'      4. cumulative rf value
#' @export
processConcordeOutput <- function(filename, markerdata, rfmat)
{
  # IMPORTANT: Concorde considers the first node to be #0, so we will need to adjust the indices to match those of the input file.

  # Read in Concorde's solution file.
  path <- scan(file=filename)

  # Remove first element of list, since this is just the number of nodes.
  path = path[-1]

  # Increment all values in path by 1 because Concorde starts indexing at 0 but R starts at 1.
  path <- path + 1

  # We want to remove the dummy marker from path.  This marker always has the largest index value, which is 1 greater than the length of markerdata.
  path = path[-which(path == (length(markerdata)+1))]

  # Create a list of marker names.
  markernames = array()

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
