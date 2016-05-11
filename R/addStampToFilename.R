#' Add datestamp to file
#'
#'    This function adds current time and date to filename and returns value
#' @usage addStampToFilename(filename, extension)
#' @param filename file to get date from
#' @param extension file extension (ie .csv)
#' @return Returns filename with date and time added
#' @export
addStampToFilename = function(filename, extension){
  filenameWithStamp = paste(filename, "_", format(Sys.time(),"%Y%m%d_%H%M"), ".", extension, sep="")
  return (filenameWithStamp)
}
