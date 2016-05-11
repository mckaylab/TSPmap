
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
Concordepath = "~/TSPexecs/concorde"


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
filepath = "results/"

duplicates = finddups(rawdata, threshold = 99,
                      filename=paste0(datapath, filepath,
                                      paste("JM20dups","tsv", sep=".")))

# Remove duplicate markers from rawdata.
rawdata = removedups(rawdata, duplicates)
