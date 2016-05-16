#' Reorder raw marker data
#'
#' Reorders raw markerdata and adds linkage group and marker position based on TSPmap results, if written as .csv can be imported to R/qtl by qtl:read.cross()
#'  @usage reorderMarkers(rawdata, distancedf)
#'  @param rawdata markerdata used to create map
#'  @param distancedf dataframe containing markernames (markername) distances (`Haldane (cM)`) and groups (groupname) of markers
#'  @return data.frame of markers properly ordered with linkage group and position
#'  @export
reorderMarkers<-function(rawdata, distancedf){
orderedmarkers<-as.matrix(rawdata[,as.character(distancedf$markername)])
orderedmarkers<-data.frame(orderedmarkers, stringsAsFactors=FALSE)
orderedmarkers<-rbind(as.character(distancedf$groupname), as.character(distancedf$`Haldane (cM)`), orderedmarkers)
return(orderedmarkers)
}
