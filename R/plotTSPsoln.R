#' Plot the TSP solution
#'
#'    Plot the TSP solution.
#' @usage plotTSPsoln(tourData, strainName, solverName, colorvector)
#' @param tourData - TSP solution object, as produced by TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput
#' @param strainName - name of the strain, as a string
#' @param solverName - name of the solver (Concorde or LKH), as a string
#' @param colorvector (OPTIONAL) - vector of color values to be used in the plot
#' @export
plotTSPsoln <- function(tourData, strainName, solverName, colorvector=NULL)
{

  if(length(colorvector) == 0)
  {
    plot(tourData[,4], cex=.5, main = paste(solverName, " Solution - ", strainName), xlab = "Marker Position", ylab = "Cumulative Recombination Fraction")
  } else {
    par(cex.lab = 0.8, cex.main = 1.2)
    plot(tourData[,4], cex=.5, col = colorvector, main = "Hamiltonian Path - IR64B", xlab = "Marker Position", ylab = "Cumulative Recombination Fraction", pch = 19)
  }

}
