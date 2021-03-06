% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findZeros.R
\name{findZeroes}
\alias{findZeroes}
\title{Identify redundant markers in the solution object that have rf = 0 so that they can be removed.}
\usage{
findZeroes(TSPsolution, rfmat)
}
\arguments{
\item{TSPsolution}{- TSP solution object, produced by the TSPinterface.R/processConcordeOutput or TSPinterface.R/processLKHOutput functions}

\item{rfmat}{recombination frequency matrix}
}
\value{
List of marker indices that can be removed.
}
\description{
This function creates a list of markers in the solution marker list which have 0 recombination frequency with their neighbors, so that they can be removed.
   NOTE: after running this function, the redundant marker list need to be removed from the solution marker list, the rf matrix, and the raw data object.  This can be done with the removeMarkers function, the pareMatrix function, and the removedups function, respectively.
}

