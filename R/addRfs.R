#' Internal function to add rf and cumulative rf values to a solution object.
#'
#'    Internal function to add rf and cumulative rf values to a solution object.
#' @usage addRfs(TSPsolution, rfmat)
#' @param TSPsolution - TSP solution
#' @param rfmat - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#'  @return
#'    TSP solution object, containing:
#'      1. rf value between each marker and the next
#'      2. cumulative rf value
#' @export
addRfs <- function(TSPsolution, rfmat)
{
  rflist = array()

  # The list of marker indices is the 1st column of TSPsolution
  path = as.numeric(TSPsolution)

  for (i in 1:(length(path)-1) )
  {

    # Since rfmatrix has the recombination frequencies in its lower triangle, we can only call elements where the x coord is > than the y coord.
    if(path[i] > path[i+1])
    {
      rflist[i] = (rfmat[path[i], path[i+1]])
    }
    else
    {
      rflist[i] = (rfmat[path[i+1], path[i]])
    }
  }

  # The last element in the weight array needs to be the weight between the last marker and the first.
  rflist[length(path)] = (rfmat[path[length(path)], path[1]])

  # Create a cumulative running total of recombination frequency.
  cumulRF = array()
  cumulRF[1] = 0

  for (i in 2:length(path))
  {
    cumulRF[i] <- cumulRF[i-1] + rflist[i-1]
  }

  # Use cbind to conbine the two arrays into one 2D array.
  returnlist <- cbind(rflist,cumulRF)

  return(returnlist)
}
