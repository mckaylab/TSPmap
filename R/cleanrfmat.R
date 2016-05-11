#' Clean up rfmat object
#'
#' @description  function to clean up rf matrix and make it pretty
#' @usage cleanrfmat(rfmat)
#' @param rfmat - recombination frequency matrix
#' @return pretty rf matrix
#' @export
#'
cleanrfmat = function(rfmat){
# round it
#   rfmat = round(rfmat, 2)
# keep lower
rfmat[upper.tri(rfmat)] = NA
diag(rfmat) = NA
return(rfmat)
}
