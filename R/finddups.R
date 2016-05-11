#' This function finds duplicate markers in the main marker data object.
#'
#' @description  Only exact duplicates are identified, meaning that the markers must have the same value in each location, and also must be missing data in the same locations.
#'  To remove the duplicate markers, see the removedups function.
#' @usage finddups(markerdata, threshold, filename)
#' @param markerdata - main marker data file, produced by the readfile.R/readFile function
#' @param threshold (OPTIONAL) - similarity threshold to use to identify how similar two markers must be in order to be considered "identical".  If this is omitted, the default behavior is to find only markers which are exact duplicates, including that they must be missing data in the same spots.  Any other value omits any spot where data are missing (from either marker), so threshold = 100 will find markers that have the same calls in every spot where both markers contain data, threadhol = 99 will find markers that are 99percent similar in every spot where both markers contain data, etc.
#' @param filename (OPTIONAL) - path and output filename as a string, if the user wishes to save the list of duplicate markers to a file
#' @return list of duplicate markers
#' @export
#'
finddups <- function(markerdata, threshold = -1, filename=NULL)
{
  # Number of markers in markerdata.
  numcols = dim(markerdata)[2]

  # Convert rawdata into a numeric matrix.
  newmatrix = data.matrix(markerdata)

  # Convert all values corresponding to '-' to a zero, so we can count them as missing in the C script.  NOTE: at this point, we treat "h" calls the same as missing calls.  This may change in the future.

  # We have to do the same for values of A and B, since we can't be sure that they will always have the same level index across all markers.
  for(i in 1:(numcols))
  {
    newmatrix[,i] = replace(newmatrix[,i], which(markerdata[,i]=="-"), as.integer(0))
    newmatrix[,i] = replace(newmatrix[,i], which(markerdata[,i]=="a"), as.integer(1))
    newmatrix[,i] = replace(newmatrix[,i], which(markerdata[,i]=="b"), as.integer(2))
    newmatrix[,i] = replace(newmatrix[,i], which(markerdata[,i]=="h"), as.integer(0))
  }

  dups = .Call("findDups", newmatrix, threshold, PACKAGE = "TSPmap")

  # The findDups() function passes back a matrix which holds the value of the original marker in each element.  In other words, if marker 120 is found to be a duplicate of marker 15, then cell 120 will have the value of 15.
  # This is so that when we find duplicates, we can know which markers they are duplicates of.
  # Therefore, the indices that we need to remove from markerdata are the indices of dups which are not -1, but not the values stored in those indices.

  # To find the indices which are not -1, we use which(dups != -1).

  # The dups list has -1 values in each element that is not a duplicate.
  # Use the unique() function because there will be many values of -1.
  # We sort the resultant list in decreasing order, which will make it easier to remove the elements from the marker list.
  duplist = sort(unique(which(dups != -1)[1:length(dups)]), decreasing = TRUE)

  writeLines(paste('\n\nfound',length(duplist),'duplicate markers'))

  # Create a list of duplicate marker names.
  # First create a copy of duplist, but sorted in increasing order, so the list of marker names will be in the correct order.
  droppedMarkers = sort(duplist)
  # Select the colnames that correspond to these indices.
  droppedMarkers = colnames(markerdata)[droppedMarkers]

  origMarker = vector()
  dupMarker = vector()

  for(i in 1:length(dups))
  {
    if(dups[i] != -1)
    {
      #print(i)
      #print(dups[i])
      origMarker = c(origMarker, colnames(markerdata)[dups[i]])
      dupMarker = c(dupMarker, colnames(markerdata)[i])
    }
  }

  duplicateMarkers = data.frame(origMarker, dupMarker, stringsAsFactors=TRUE)

  # Sort the duplicateMarkers list by the first column and write it to a file.
  newdups = duplicateMarkers[with(duplicateMarkers, order(origMarker)),]

  # Write this out to a file.
  if(!is.null(filename))
    write.table(newdups, file = filename, row.names=FALSE, sep="\t")

  return(duplist)
}
