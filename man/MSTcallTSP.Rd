% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MSTcallTSP.R
\name{MSTcallTSP}
\alias{MSTcallTSP}
\title{Call the LKH solver and process the results.}
\usage{
MSTcallTSP(MSTgroup, rawdata, rfmat, LKHexec)
}
\arguments{
\item{MSTgroup}{cluster of marker indices.}

\item{rawdata}{raw marker data.}

\item{rfmat}{rf matrix.}

\item{LKHexec}{path and filename of LKH executable.  MUST NOT contain the tilde character.}
}
\value{
List of marker IDs in order of TSP solution.
}
\description{
This function calls the LKH solver.
}

