
# function to clean up rf matrix and make it pretty
cleanrfmat = function(rfmat){
  # round it
#   rfmat = round(rfmat, 2)
  # keep lower
  rfmat[upper.tri(rfmat)] = NA
  diag(rfmat) = NA
  return(rfmat)
}

# function to add marker names
addMarkerNames = function(mymatrix, rawdata){
  # add row names after making df
  mydf = as.data.frame(mymatrix)
  # add colnames
  colnames(mydf) = names(rawdata)
  mydf$marker = names(rawdata)

  # swap col order
  mydf = cbind(mydf[,ncol(mydf)], mydf[,1:ncol(mydf) - 1])
  colnames(mydf)[colnames(mydf) == colnames(mydf)[1]] = "marker"
  return(mydf)
}
