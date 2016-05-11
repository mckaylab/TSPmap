#' Create TSP problem file.
#'
#'  Create TSP problem file (for use with both Concorde and LKH).
#'  This converts the rf matrix to an integer matrix that the TSP solvers need.
#'  HAMILTONIAN PATH CORRECTION: This also adds a dummy node which has 0 weight to all other nodes so that the TSP solver will find a Hamiltonian path instead of a circuit.
#' @usage createTSPFile(markerdata,rfmatrix,strainName, dirname)
#' @param markerdata main marker data file, produced by the readfile.R/readFile function
#' @param rfmatrix recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param strainName name of the strain as a string (no spaces), used to name the file
#' @param dirname directory in which to create the file
#' @return Returns the path and filename of the TSP file (.tsp extension).
#' @export
createTSPFile <- function(markerdata, rfmatrix, strainName, dirname)
{
  # Disable scientific notation so the .tsp files does not have character strings in it.
  options(scipen=999)

  # Update numcols with the new value, since we have removed duplicate markers.
  numcols = dim(markerdata)[2]

  # All values need to be turned into ints so that LKH can read them in, so multiply by 100000 so that we can preserve 5 significant digits (i.e. 0.12345 becomes 12345).  Any factor larger than this seems to cause problems in LKH.
  tspmatrix <- round(rfmatrix*100000)

  # Create data file that will be sent to Concorde.
  filename = paste0(strainName, ".tsp")

  # Add dirname.
  filename = paste0(dirname, filename)

  # First output the header.
  # append=FALSE by default, so this first line clears out any existing file data.
  write(paste("NAME : ", strainName, sep=""), file = filename)

  write("TYPE : TSP", file = filename, append = TRUE)
  write("COMMENT : Recombination frequency matrix of bio markers.", file = filename, append = TRUE)
  # Add 1 to the size of the matrix because of the dummy marker that we will add later.
  write(paste("DIMENSION :",dim(rfmatrix)[1]+1, collapse=""), file = filename, append = TRUE)
  write("EDGE_WEIGHT_TYPE: EXPLICIT", file = filename, append = TRUE)
  write("EDGE_WEIGHT_FORMAT: UPPER_ROW", file = filename, append = TRUE)
  write("EDGE_WEIGHT_SECTION", file = filename, append = TRUE)

  # End at dim(tspmatrix)[1]-1 because there's nothing to output for the last row (since we are not outputting the diagonal elements).

  for(i in 1:(dim(tspmatrix)[1]-1))
  {
    # Start at i+1 to exclude the diagonal elements.
    # Also append a zero at the end to represent the dummy node.
    outputthis = c(tspmatrix[(i+1):dim(tspmatrix)[1],i],0)
    write(outputthis, file = filename, append = TRUE, ncolumns=length(outputthis))
  }

  # Apply Hamilton path modification.
  # Instead of having the TSP solver find a Hamiltonian circuit, we will add a dummy node with rf = 0 at all elements so that the TSP solver will produce a path that does not have to connect head-to-tail.
  # NOTE: this marker is removed from the TSP solution in the processLKHOutput/processConcordeOutput methods.

  # Add a final zero for the dummy node.
  write(0, file = filename, append = TRUE, ncolumns=1)
  write("EOF", file = filename, append = TRUE)

  # Reset scientific notation to the default value.
  options(scipen=0)

  return(filename)
}
