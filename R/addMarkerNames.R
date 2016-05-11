#' add marker names
#'
#' @description  function to add marker names
#' @usage finddups(markerdata, threshold, filename)
#' @param mymatrix - rfmatrix
#' @param rawdata raw marker data file to get names from
#' @return pretty rf matrix
#' @export
#'
addMarkerNames = function(mymatrix, rawdata){
  # add row names after making df
  mydf = as.data.frame(mymatrix)
  # add colnames
  colnames(mydf) = names(rawdata)
  mydf$marker = names(rawdata)

  # swap col order
  mydf = cbind(mydf[,ncol(mydf)], mydf[,1:ncol(mydf) - 1])
  colnames(mydf)[colnames(mydf) == colnames(mydf)[1]] = "marker"
  return(mydf)
}
