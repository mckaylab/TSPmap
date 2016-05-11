#' Cap the rf matrix according to the given threshold.
#'
#'    For all recombination frequency value that are above the threshold parameter, this function raises these values to capval.  This is done in order to reduce the runtime of the TSP solvers.
#' @usage capRFmat(rfmat, threshold, capval)
#' @param rfmat - recombination frequency matrix.
#' @param threshold - upper limit on meaningful rf values.
#' @param capval - value to which rf values above threshold will be set.
#' @return matrix of recombination frequency values, with all values above threshold set to capval
#'
capRFmat <- function(rfmat, threshold, capval = 1.0)
{
  result = ifelse(rfmat > threshold, capval, rfmat)

  return(result)
}
