#' Write the sorted RF matrix to file.
#'
#' This function saves a sorted RF matrix to file, along with the marker names.
#' @usage writeSortedRFmat(filename, pathData, sortedrfmatrix)
#' @param filename - filename to create
#' @param pathData - TSP solution with marker names
#' @param sortedrfmatrix - rf matrix in sorted order
#' @return matrix of recombination frequency values with marker names added
#' @export
writeSortedRFmat <- function(filename, pathData, sortedrfmatrix)
{
  # Want to add marker names to first column and row to result here so it's easier to look at in Excel.
  rfmatrix = cbind(pathData[,2], sortedrfmatrix)
  rfmatrix = rbind(c(0,pathData[,2]), rfmatrix)

  # Used to use write.csv here, but write.table gives us more control.
  write.table(rfmatrix, file = filename, sep=",", col.names=FALSE, row.names=FALSE)

  return(rfmatrix)
}
