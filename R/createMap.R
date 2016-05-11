#' Process a linkage group with Concorde to get the final linkage map.
#'
#'  Starting with the clusters produced by autoClusterMST, this function merges them until we are left with full linkage groups.
#' @usage createMap(rawdata, rfmat, cluster, Concexec)
#' @param rawdata - raw marker data
#' @param rfmat - recombination frequency matrix, produced by the rfmatrix.R/computeRFmat function
#' @param cluster - one of the linkage groups produced by autoMergeClusters
#' @param Concexec - path and filename of Concorde executable - must not contain the tilde character!
#' @return Returns the final linkage map.
#' @export
createMap <- function(rawdata, rfmat, cluster, Concexec)
{
  # Find the submatrix of the rf matrix that includes all these markers.
  submatrix = pareMatrix(rfmat, unlist(cluster))
  subrawdata = rawdata[unlist(cluster)]

  # Create TSP file in R's temporary directory, send it to Concorde.
  # NOTE: tempdir() needs to have an extra "/" added to it.
  # Also, tempdir() appears to have an extra "/" in it, remove this before continuing.
  tempdirectory = paste0(tempdir(),"/")
  tempdirectory = gsub("//", "/", tempdirectory)

  filename = "tsptempfile"
  tspfile = createTSPFile(subrawdata, submatrix, filename, tempdirectory)

  # Run Concorde.
  solnfile_Conc = callConcorde(Concexec, tspfile, tempdirectory)

  # Process Concorde solution file.
  ConcordeMarkerlist = processConcordeOutput(solnfile_Conc, subrawdata, submatrix)

  # Now the marker indices need to be translated back to their proper values (i.e. the values within rawdata).
  for(j in 1:length(ConcordeMarkerlist[,1]))
  {
    ConcordeMarkerlist[j,1] = which(colnames(rawdata) == ConcordeMarkerlist[j,2])
  }


  return(ConcordeMarkerlist)
}

