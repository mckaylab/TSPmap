#' create table of marker distances and groups
#'
#' @description converts distancedf into table that can be converted into a map object by qtl::table2map
#' @usage distancedf2table(distancedf)
#' @param distancedf dataframe containing markernames (markername) distances (`Haldane (cM)`) and groups (groupname) of markers
#' @return table with chr=linkage group, pos=distance, and row.names = marker names
#' @export
distancedf2table<-function(distancedf){
markertable<-data.frame(chr=distancedf$groupname, pos=as.character(distancedf$`Haldane (cM)`))
rownames(markertable)<-distancedf$markername
return(markertable)
}
