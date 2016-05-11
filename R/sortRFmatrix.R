#' This function sorts the rfmatrix so that the markers are in the same order as the TSP tour.
#'
#'    This function sorts the rfmatrix so that the markers are in the same order as the TSP tour.
#'    The resultant matrix of this function is used in the cutpoints.R/matchGroups function.
#'    NOTE: The resultant matrix of this function cannot be sent to the writeRFmat function because the columns/rows no longer match up with the main markerdata object.
#' @usage sortRFmatrix(rfmatrix, pathData, saveFlag)
#' @param rfmatrix - matrix of recombination frequency values produced by the rfmatrix.R/computeRFmat function
#' @param pathData - TSP tour and weights, produced by the TSPinterface.R/processConcordeOutput function
#' @param saveFlag - if saveFlag == 1, save sorted matrix to file
#' @return matrix of recombination frequency values sorted in order of TSP tour
#' @export
sortRFmatrix <- function(rfmatrix, pathData, saveFlag=0)
{
  # Pull out the first column of the pathData object and convert into integers.
  path = as.integer(pathData[,1])

  # Make a copy of result so we don't mess with the original.
  sortedRF = rfmatrix

  # For tracking purposes, replace the diagonal elements with the marker index (which to start off will just be in order).  When we are done, the diagonals should be in the same order as shown in path.

  for(i in 1:length(path))
    sortedRF[i,i] = -i

  # Keep track of where each marker's weights are, since we will be moving them around.
  # Each ith element of this list holds the row/col index where the ith marker is currently being held.
  # eg, location[15] = 6 means that marker #15 is in row/col #6.
  location = 1:length(path)

  # As we step through the rows/cols of the sortedRF matrix, the ith element of path shows us which marker should be in position i, and the [path]th element of location shows us where that marker is currently located.  In short, we want to move the sortedRF[location[path[i]]] to sortedRF[i].

  # Also, row/col #1 of the result matrix are just zeroes so row/col #2 correspond to marker #1.  Because of this, we need to adjust the indices by 1 when referencing sortedRF.

  limit = length(path)
  for (i in 1:limit)
  {
    # Initially, the rf matrix holds the markers in order of their index value, that is, marker #1 is in row/col #1, marker #2 is in rol/col #2, etc.
    # As we loop over i, we want to put the ith element of path into row/col i.
    # eg, for i = 6, we find which marker is at path[6] and move it from wherever it's currently stored in the rf matrix into row/col 6.
    # Since things get swapped around, we can't just assume that row/col 6 holds marker #6.  Instead we have to find which element of location[] has the value 6, and the corresponding index will be the value of the marker which is stored there.
    # Call this moveMarker, because this marker is just getting moved around.
    moveMarker = which(location == i)

    # Now we swap row/col i with the row/col where marker path[i] is currently stored.
    # Call this fixMarker, because once we are done, this marker will be fixed in its correct spot.

    fixMarker = location[path[i]]

    # If we have removed the rf = 0 markers, we could get some NA values here.  Only proceed if we do not have a NA value.
    if(!is.na(fixMarker))
    {
      # First swap row i with row fixMarker.
      temp = sortedRF[fixMarker,]
      sortedRF[fixMarker,] = sortedRF[i,]
      sortedRF[i,] = temp

      # Now swap the columns.
      temp = sortedRF[,fixMarker]
      sortedRF[,fixMarker] = sortedRF[,i]
      sortedRF[,i] = temp

      # At this point, we have put marker #fixMarker into position i, and marker #moveMarker into position fixMarker.
      location[path[i]] = i

      location[moveMarker] = fixMarker
    }
  }

  # Write out the matrix to a file.
  if(saveFlag == 1)
  {
    print("Saving sorted RF matrix")
    write.csv(sortedRF, file = "sortedRFmatrix.csv", row.names=FALSE)
  }

  return(sortedRF)
}


