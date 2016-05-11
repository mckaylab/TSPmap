
# Get TSP solvers ---------------------------------------------------------

# http://www.akira.ruc.dk/~keld/research/LKH/
#   http://www.akira.ruc.dk/~keld/research/LKH/LKH-2.0.7.tgz
# tar xvfz LKH-2.0.7.tgz
# cd LKH-2.0.7
# make
#
# http://www.math.uwaterloo.ca/tsp/concorde/downloads/downloads.htm
# download and unzip executable file.
# chmod u+x concorde
#http://www.math.uwaterloo.ca/tsp/concorde/DOC/README.html

# Install TSPmap ----------------------------------------------------------

#install.packages(devtools)
library(devtools)
install_bitbucket("greymonroe/TSPmap")
library(TSPmap)

# Setting up --------------------------------------------------------------

# Define your path to the data directory
# this is where input and output files are located
# datapath <- "/path/to/data/"
datapath <- "~/Desktop/data/"

# within this directory, we recommend the following directories:
# /rawdata/
# /results/LKH
# /results/Concorde
# /TSPfiles/

# Define poath ot LKH and Concorde exec files.
LKHexec = "~/TSPexecs/LKH"
Concordepath = "~/TSPexecs/concorde_mac"

# -----------------------------------------------------
# 1 read in data file
# -----------------------------------------------------

# format can be either comma separated or tab delimited
# The data MUST contain ONLY the characters "a", "b", "h", and "-", to indicate
# missing calls.
# rows are individuals in the population, columns are markers
# If the second column isn't a "Classification" column for joinmap,
# then set the flag to false:
rawdata = readFile(paste0(datapath, "rawdata/JM20Demo_rawinput.tsv"),
                   classificationFlag=F)

# -----------------------------------------------------
# 2 find and remove duplicate markers
# -----------------------------------------------------

# If you don’t want to save the duplicate list to a file, just do this:
duplicates = finddups(rawdata, threshold = 99)

# If you also want a list of duplicates exported, you can do this instead:

# first, add timestamp to file
# don't do this for the JM20 data
# filename = addStampToFilename("JM20Dups")
# define where you want it saved
resultpath = "results/"

duplicates = finddups(rawdata, threshold = 99,
                      filename=paste0(datapath, resultpath,
                                      paste("JM20dups","tsv", sep=".")))

# Remove duplicate markers from rawdata.
rawdata = removedups(rawdata, duplicates)

# -----------------------------------------------------
# 3 calculate recombination frequency matrix
# -----------------------------------------------------

# Compute recombination frequency matrix.

# We will use the optional correction flag parameter.  This flag slightly
# modifies the recombination frequency matrix to avoid situations where markers
# with missing calls might act as a bridge between other markers, resulting in
# false correlations between markers which should not be related.
# This does not affect runtime but may affect solution quality, especially in
# data sets with many missing calls.

rfmat = computeRFmat(rawdata, corrFlag=T)

# it is easier to look at if you clean it up:
rfmatclean = cleanrfmat(rfmat)
rfmatclean = addMarkerNames(rfmatclean, rawdata)

# OPTIONAL: Write out the RF matrix to a file after running computeRFmat.
# define what you want to call the file
# filename = addStampToFilename("JM20RFmatrixFromTSPmap")
filepath = "results/"

# you can't use the rfmatclean here, though obviously you could write rfmatclean to a file as well
writeRFmat(paste0(datapath, filepath, paste("JM20RFmatrixFromTSPmap","csv", sep=".")), rawdata, rfmat)

# Now we set an upper limit on the meaningful rf values.  Any recombination
# frequency values above this threshold are too large to contribute to the
# formation of the linkage groups, so they are effectively equivalent and will
# be replaced with a value of 0.5.  This allows the TSP solver to run more
# efficiently in order to reduce runtime.
# If this threshold is set too low, useful information may be lost and the
# linkage groups may be incorrect.  If this threshold is set too high, the TSP
# solver will take longer to run.

rfCap = 0.4
cappedrfmat = capRFmat(rfmat, rfCap, 0.5)

# The capped RF matrix will be used in all computations from this point onward.

# -----------------------------------------------------
# 4 linkage group formation
# -----------------------------------------------------

# define how many chromosomes in your species
numChromosomes = 5

# Break into small clusters.  This is the first step in forming linkage groups and ensures that the markers from different linkage groups are all well-separated.

# The internalRfThreshold variable sets the upper limit on the recombination
# frequency values that will be allowed within a cluster of markers.  If this value is set too high, the linkage groups may not be well-separated and the linkage groups will be incorrect.  If this value is set too low, many small clusters will be formed resulting in longer run times.

rawclusters = autoClusterMST(rawdata, cappedrfmat, numChromosomes, LKHexec, internalRfThreshold=0.4)

# warning is ok:
# Warning message:
#   In clusterOrder[i] = which(clusterSizes == sortedSizes[i]) :
#   number of items to replace is not a multiple of replacement length

# then merge the clusters into the final linkage groups.  This will create a number
# of clusters equal to the value of numChromosomes.

mergedclusters = autoMergeClusters(rawclusters, cappedrfmat, LKHexec, numChromosomes, rfCap)

# for each linkage group, compute the final linkage map.

finalMap = list()
for(i in 1:length(mergedclusters))
{
  finalMap[[i]] = createMap(rawdata, cappedrfmat, mergedclusters[i], Concordepath)
}

# remove empty elements, if necessary
finalMap = finalMap[!sapply(finalMap, is.null)]


# give the groups better names (it gives them all "1" by default. should prob fix this)
# NOT SURE WHAT THIS MEANS??? -ZAA

# also within this loop it calculates the genetic distance. can do kosambi instead
newdf = data.frame()
for(i in 1:length(finalMap)){
  # print (i)
  temp = as.data.frame(computeHaldane(finalMap[[i]]))
  temp$groupname = i
  newdf = rbind(newdf, temp)
}

# -----------------------------------------------------
# 5 examine whether these groups make sense
# -----------------------------------------------------

# if a group is in the wrong order, you can reverse it:
# here, we are reversing group 4
finalMap[[4]] = reverseCluster(finalMap[[4]], rfmat)

# Define a list of colors to be used on plots, if desired.
colors = c("red1", "green", "lightblue", "blue", "orange")

# These variables will be used to create an overall plot of all linkage groups.
fullMap = list()
fullColorVec = list()

# Create individual plots of the cumulative rf for each linkage group.
for(i in 1:length(finalMap))
{
  # Create a color vector for each plot:
  colorVec = rep(colors[i], length(finalMap[[i]][,1]))

  plotTSPsoln(finalMap[[i]], "JM20", "Concorde", colorVec)

  # This saves the plots.
  filename = paste0(datapath, resultpath, "/JM20",i,".png")
  dev.copy(png,filename,width=8,height=8,units="in",res=100)
  dev.off()

  # Update the lists for the overall plot.
  fullMap = c(fullMap, as.numeric(finalMap[[i]][,4]))
  fullColorVec = c(fullColorVec, list(colorVec))
}

fullMap = unlist(fullMap)
fullColorVec = unlist(fullColorVec)

# Create the overall plot.
plot(fullMap, cex=.5, col = fullColorVec, xlab = "Marker Position", ylab = "Marker Position in JoinMap Solution", pch = 19)




