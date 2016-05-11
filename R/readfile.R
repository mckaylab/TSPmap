#' This function reads in a raw data file from disk.
#'
#' This function reads in a raw data file from disk. This is the only place where data are brought into the package, all calculations are done from these data.
#' File format:
#'    1. Tab-separated or comma-separated file
#'    2. Markers are arranged in rows
#'    3. One header row which can be discarded
#'    4. One header column with the marker names as strings
#'    5. OPTIONAL: A second column which can also be discarded (this column typically has the string "(a,b)" in each cell).  If this column is not present, set the classificationFlag parameter to FALSE.
#' @usage readFile(filename,classificationFlag,transpose)
#' @param filename - path and input filename as a string
#' @param classificationFlag - flag to indicate whether the 2nd column in the data file is the classification column, which will be discarded.  Default is TRUE.
#' @param transpose - flag to indicate whether or not the data needs to be transposed.  Default is TRUE.
#' @return main marker data object
#' @export
readFile <- function(filename, classificationFlag=TRUE, transpose = TRUE)
{
  # Check the file format (tsv or csv) and import the file appropriately.
  filenamelength = nchar(filename)

  if(substr(filename, (filenamelength-2), filenamelength) == "csv")
    data <- read.table(file=filename, sep=",", na.strings = "")
  else
    data <- read.table(file=filename, sep="\t", na.strings = "")

  # If classificationFlag == true, column #2 is useless so get rid of it.
  if(classificationFlag)
    data[[2]] <- NULL

  if(transpose)
  {
    # Take the transpose.
    data = t(data)

    # Now the first column is not needed.
    data <- data[,-1]

    # Set column names equal to the first column, then delete the first column.
    colnames(data) = data[1,]

    # Get rid of the first row.
    data <- data[-1,]
  } else {
    # Set column names equal to the first row, then delete the first row.
    colnames(data) = sapply(data[1,], as.character)
    data <- data[-1,]
  }

  # Row names are not needed.
  rownames(data) = NULL

  # Need to make a data frame out of this, so that the values are actually stored numerically, because this is what finddups is expecting.
  data = as.data.frame(data)

  return(data)
}
