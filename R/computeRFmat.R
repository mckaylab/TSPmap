#' This function computes the recombination frequency matrix.
#'
#'    This function computes the recombination frequency matrix.
#' @usage computeRFmat(markerdata, corrFlag)
#' @param markerdata - main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function.
#' @param corrFlag - flag which controls whether or not the missing data correction method is used.  Default is FALSE.
#' @return matrix of recombination frequency values
#' @export
computeRFmat <- function(markerdata, corrFlag = FALSE)
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

  if(corrFlag)
  {
    datacorrection = 1
  } else {
    datacorrection = 0
  }

  # Call the C code to compute the rf matrix, and round the result to 5 decimal places.  The reason we do this is because it seems that the weight matrix in LKH cannot handle weight values larger than 5 digits, so we don't lose any information by doing this.
  # Also see TSPinterface.R/createTSPFile
  result = round(.Call("computeRF", newmatrix, datacorrection, PACKAGE = "TSPmap"), digits=6)

  # Put the number of valid individuals in the diagonal slots.
  # We need this so that when we look for similar rf vectors, we know which marker to keep (the one with more valid individuals).
  for(i in 1:(numcols))
  {
    result[i,i] = 0
  }

  return(result)
}
