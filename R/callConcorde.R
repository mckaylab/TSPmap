#' Run Concorde TSP solver.
#'
#'    Run Concorde TSP solver.
#' @usage callConcorde(execpath,filename, outputdir)
#' @param execpath - path and name of Concorde executable.
#' @param filename path and input filename as a string
#' @param outputdir - directory in which Concorde will save the solution file
#' @return Returns the path and filename of the Concorde solution file (.sol extension).
#' @export
callConcorde <- function(execpath, filename, outputdir)
{
  # Strip out the base filename.
  locs = gregexpr("/", filename)
  firstspot = locs[[1]][length(locs[[1]])]+1

  firstsubstring = substring(filename, firstspot)

  locs2 = gregexpr("\\.", firstsubstring)
  secondspot = locs2[[1]][1] - 1

  fileroot = substring(firstsubstring, 1, secondspot)

  # Replace all spaces with "\\\\ "
  execpath = gsub(" ", "\\\\ ", execpath)
  filename = gsub(" ", "\\\\ ", filename)
  outputdir = gsub(" ", "\\\\ ", outputdir)

  # Create the output file name (.sol file extension).
  outfilename = paste0(outputdir, fileroot, ".sol")

  # Need to add these parameters to the concorde call:
  #   -x    delete files on completion (sav pul mas)
  #   -s n  seed value
  #   -o f  output file name (for optimal tour)
  exec = paste(execpath, '-x -s 0 -o', outfilename)

  system(paste(exec, file=filename), wait=TRUE)

  return(outfilename)
}
