#' Call the LKH solver and process the results.
#'
#' This function calls the LKH solver.
#' @usage MSTcallTSP(MSTgroup, rawdata, rfmat, LKHexec)
#' @param MSTgroup cluster of marker indices.
#' @param rawdata raw marker data.
#' @param rfmat rf matrix.
#' @param LKHexec path and filename of LKH executable.  MUST NOT contain the tilde character.
#' @return List of marker IDs in order of TSP solution.
#' @export
MSTcallTSP = function(MSTgroup, rawdata, rfmat, LKHexec)
{
  cluster = as.numeric(unlist(MSTgroup))

  # Find the submatrix of the rf matrix that includes all these markers.
  submatrix = pareMatrix(rfmat, cluster)
  subrawdata = rawdata[cluster]

  # Create TSP matrix, send it to LKH.
  # NOTE: tempdir() needs to have an extra "/" added to it.
  tempdirectory = paste0(getwd(),"/")
  filename = "tsptempfile"
  tspfile = createTSPFile(subrawdata, submatrix, filename, tempdirectory)

  solnfile_LKH = callLKH(LKHexec, tspfile, tempdirectory)

  # Process the results, see what the largest weight is.
  SubMarkerlist = processLKHOutput(paste0(tempdirectory,filename,".LKH"), subrawdata, submatrix)

  fixedmarkerlist = cluster[as.numeric(SubMarkerlist[,1])]

  SubMarkerlist[,1] = fixedmarkerlist

  return(SubMarkerlist)
}
