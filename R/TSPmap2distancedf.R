#' convert final tspmap to distance data.frame
#'
#' converts final map created by concorde to a dataframe including haldane distances and groupname
#' @usage TSPmap2distancedf(finalMap)
#' @param finalMap final map list object created with TSPmap after grouping and ordering markers
#' @return dataframe containing columns: index, markername, rf, cumulative rf, Haldane (cM), group, groupname
#' @export
#'
TSPmap2distancedf<-function(finalMap){
distancedf = data.frame()
for(i in 1:length(finalMap)){
  # print (i)
  temp = as.data.frame(computeHaldane(finalMap[[i]]))
  temp$groupname = i
  distancedf = rbind(distancedf, temp)
}
return(distancedf)
}
