% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkNongroupEdges.R
\name{checkNongroupEdges}
\alias{checkNongroupEdges}
\title{Find the smallest edge weight between the target cluster and all other clusters.}
\usage{
checkNongroupEdges(rfmatrix, targetIndex, groupList)
}
\arguments{
\item{rfmatrix}{- rf matrix}

\item{targetIndex}{- index of the target cluster}

\item{groupList}{- list of all clusters}
}
\value{
List of minimum edge weights between the cluster at targetIndex and all other clusters.
}
\description{
This function returns a list of the smallest edge weight between a given cluster and all other clusters.  This is used to make the decision to merge clusters when the smallest edge weight is below a given threshold.
}

