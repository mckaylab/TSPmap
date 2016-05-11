#' This function removes duplicate markers from the main marker data object.
#'
#'    This function removes duplicate markers from the main marker data object.
#' @usage removedups(markerdata, duplist)
#' @param markerdata - main marker data file, produced by the readfile.R/readFile function.
#' @param duplist - list of duplicate markers, produced by the duplicatemarkers.R/finddups function.
#' @return main marker data file with duplicate markers removed
#' @export
removedups <- function(markerdata, duplist)
{
  for (i in duplist)
  {
    markerdata[[i]] <- NULL
  }

  writeLines(paste('\n\nremoved',length(duplist),'duplicate markers, leaving', length(markerdata[1,])))

  return(markerdata)
}
