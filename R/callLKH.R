#' Run LKH TSP solver.
#'
#'    Run LKH TSP solver.
#' @usage callLKH(execpath,filename,outputdir)
#' @param execpath - path and name of LKH executable
#' @param filename - path and input filename as a string
#' @param outputdir - directory in which LKH will save the solution file
#' @return
#'    Returns the path and filename of the LKH solution file (.LKH extension).
#' @export
callLKH <- function(execpath, filename, outputdir)
{
  # Get home directory.
  home = system("echo ${HOME}", intern =TRUE)

  # Replace any ~ in the parameters with home.
  execpath = gsub("~", home, execpath)
  filename = gsub("~", home, filename)
  outputdir = gsub("~", home, outputdir)

  # Need to create LKH input file (.par), which references the TSP file created by the call to createTSPFile.
  #'
  # Strip out the base filename.
  locs = gregexpr("/", filename)
  firstspot = locs[[1]][length(locs[[1]])]+1

  firstsubstring = substring(filename, firstspot)

  locs2 = gregexpr("\\.", firstsubstring)
  secondspot = locs2[[1]][1] - 1

  fileroot = substring(firstsubstring, 1, secondspot)

  # Create the LKH input file (.par file extension).
  parfilename = paste(outputdir, fileroot, ".par", sep="")

  # First output the header.
  # append=FALSE by default, so this first line clears out any existing file data.
  write(paste("PROBLEM_FILE = ", filename, sep=""), file = parfilename)

  write("RUNS = 1", file = parfilename, append = TRUE)

  LKHoutputfilename = paste(outputdir, fileroot, ".LKH", sep="")
  write(paste("OUTPUT_TOUR_FILE = ", LKHoutputfilename, sep = ""), file = parfilename, append = TRUE)

  # If a file already exists with parfilename, delete it before running the TSP solver.  This prevents the user from thinking that the solver completed successfully in the case where the solver actually fails.
  if(file.exists(LKHoutputfilename))
    file.remove(LKHoutputfilename)

  # Replace all spaces with "\\\\ " before performing the system calls.
  execpath = gsub(" ", "\\\\ ", execpath)
  filename = gsub(" ", "\\\\ ", filename)
  outputdir = gsub(" ", "\\\\ ", outputdir)
  parfilename = gsub(" ", "\\\\ ", parfilename)

  system(paste(execpath, file=parfilename), wait=TRUE)
  return(LKHoutputfilename)
}
