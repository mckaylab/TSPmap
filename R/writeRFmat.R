#' writes the recombination frequency matrix to a file
#'
#'    This function writes the recombination frequency matrix to a file, along with the marker names in both the header row and header column.
#' File format:
#'    1. Comma-separated file
#' @usage writeRFmat(filename, markerdata, rfmatrix)
#' @param filename - path and output filename as a string
#' @param markerdata - main marker data file, produced by either the readfile.R/readFile function or the duplicatemarkers.R/removedups function
#' @param rfmatrix - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @return Rf matrix with marker names in the header row/column.
#' @export
writeRFmat <- function(filename, markerdata, rfmatrix)
{
  # Want to add marker names to first column and row to result here so it's easier to look at in Excel.
  rfmatrix = cbind(colnames(markerdata), rfmatrix)
  rfmatrix = rbind(c(0,colnames(markerdata)), rfmatrix)

  # Used to use write.csv here, but write.table gives us more control.
  write.table(rfmatrix, file = filename, sep=",", col.names=FALSE, row.names=FALSE)
}
