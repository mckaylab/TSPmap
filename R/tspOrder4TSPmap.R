#' @title Order markers within linkage groups
#'
#' @description
#' \code{tspOrder} Apply a traveling salesperson problem solver to find the
#' shortest path through the recombination fraction matrix. Requires that
#' the TSP R package and the TSP solver "Concorde" program are installed. See details.
#'
#' @param cross The QTL cross object.
#' @param concorde_path Required. The directory containing concorde executables.
#' @param return If "cross" is specified, pass the marker order through newLG, which quickly
#' reorders the markers of the original cross to follow the TSP order. Otherwise,
#' return a named list of markers in order for each chromosome.
#' @param ... Additional arguments passed on to solve_TSP that are then paseed to TSP::contol.
#' @details This function relies on the TSP R packages to perform the TSP solvers. See documentation
#' therein, esspecially the function TSP::solve_TSP. We permit inference of marker order within a
#' Hamiltonian circuit by adding a dummy node that has 0 distance to all other nodes. This
#' allows for discrete start and end points and is more appropriate for genetic map construction
#' than forcing a complete route through all markers.
#'
#' To install the concorde program:
#' http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz.
#' If on a mac, this is not super easy.
#' See https://qmha.wordpress.com/2015/08/20/installing-concorde-on-mac-os-x/
#' for details on the best way to do this.
#'
#' @return Either a cross object with reordered markers, or a named list
#' of markers in a new order.
#'
#' @examples
#' \dontrun{
#' library(qtlTools)
#' set.seed(42)
#' map<-sim.map(len = 200, n.mar = 200, include.x = F, eq.spacing=T)
#' cross<-sim.cross(map, type = "riself", map.function = "kosambi", error.prob = 0)
#' cross<-est.rf(cross)
#' markerlist<-lapply(chrnames(cross), function(x) sample(markernames(cross, chr =x)))
#' names(markerlist)<-as.character(chrnames(cross))
#' cross.rand<-newLG(cross, markerList = markerlist, keep.rf=T)
#' set.seed(42)
#' cross.ord<-tspOrder2(cross = cross.rand,
#'   concorde_path = "/Users/John/Documents/concorde/TSP")
#'
#' plot(match(markernames(cross), markernames(cross.ord)),
#' xlab = "position in similated map", ylab = "position in ordered map")
#' }
#' @import qtl
#' @import TSP
#' @export
tspOrder2<-function(cross,
                   concorde_path,
                   return = "cross"){
  if(!requireNamespace("TSP", quietly = TRUE)){
    stop("install the TSP package to use tspOrder\n")
  }else{
    requireNamespace("TSP", quietly = TRUE)
  }
  if(is.null(concorde_path)){
    stop("if method = concorde, concorde_path must specify directory of program\n")
  }

  newLG<-function(cross, markerList, keep.rf = TRUE){
    if(any(is.null(names(markerList)))){
      names(markerList)<-as.character(1:length(markerList))
    }
    # drop markers not in markerList
    newmars<-unlist(markerList)
    if(any(duplicated(newmars))){
      stop("duplicated markers found in markerList, all markers must be unique\n")
    }

    if (!("rf" %in% names(cross)) & keep.rf) {
      warning("Running est.rf.")
      cross <- est.rf(cross)
    }

    if(any(!newmars %in% markernames(cross))){
      stop("some markers in list are not in the cross, dropthem\n")
    }

    if(any(!markernames(cross) %in% newmars)){
      todrop<-markernames(cross)[!markernames(cross) %in% newmars]
      cross<-drop.markers(cross, markers = todrop)
    }

    n.mar <- nmar(cross)
    tot.mar <- totmar(cross)
    if(keep.rf) {
      rf <- cross$rf
      diagrf <- diag(rf)
      if (ncol(rf) != tot.mar)
        stop("dimension of recombination fractions inconsistent with no. markers in cross.")
      onlylod <- attr(cross$rf, "onlylod")

      lod <- rf
      lod[lower.tri(rf)] <- t(rf)[lower.tri(rf)]
      rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
      diag(rf) <- 1
      diag(lod) <- 0
    }

    marnam <- markernames(cross)
    chrstart <- rep(names(cross$geno), n.mar)
    ingrp <- 1:tot.mar
    chrnum<-1:length(markerList)
    revgrp <- rep(chrnum,sapply(markerList, length))

    cross <- clean(cross)
    chrtype <- rep(sapply(cross$geno, class), n.mar)
    crosstype <- class(cross)[1]
    g <- pull.geno(cross)
    cross$geno <- vector("list", max(revgrp))
    names(cross$geno) <- 1:max(revgrp)

    for (i in 1:max(revgrp)) {
      cross$geno[[i]]$data <- g[, markerList[[i]], drop = FALSE]
      cross$geno[[i]]$map <- seq(0, by = 10, length = length(markerList[[i]]))
      if (crosstype == "4way") {
        cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
                                     cross$geno[[i]]$map)
        colnames(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      }else{
        names(cross$geno[[i]]$map) <- colnames(cross$geno[[i]]$data)
      }
      thechrtype <- unique(chrtype[revgrp == i])
      if (length(thechrtype) > 1){
        warning("Problem with linkage group ", i, ": A or X?\\n",
                paste(thechrtype, collapse = " "))
      }else{
        class(cross$geno[[i]]) <- thechrtype
      }
    }
    mname <- markernames(cross)
    m <- match(mname, marnam)
    if(keep.rf) {
      rf <- rf[m, m]
      lod <- lod[m, m]
      rf[upper.tri(rf)] <- lod[upper.tri(lod)]
      diag(rf) <- diagrf[m]
      cross$rf <- rf
    }
    return(cross)
  }

  concorde_path(concorde_path)

  rf<-data.matrix(pull.rf(cross, what = "rf"))
  class(rf)<-"matrix"
  diag(rf)<-0

  markerList<-lapply(chrnames(cross), function(x) markernames(cross, chr = x))
  names(markerList)<-chrnames(cross)

  chr.ord<-lapply(names(markerList),function(x){
    mnames<-markerList[[x]]
    totsp<-rf[mnames,mnames]
    chr.tsp<-TSP(totsp)

    tsp <- insert_dummy(chr.tsp, label = "cut")
    tour <- solve_TSP(tsp, method="concorde")
    path <- cut_tour(tour, "cut")

    ord<-labels(path)
    return(ord)
  })
  names(chr.ord)<-names(markerList)
  if(return != "cross"){
    return(chr.ord)
  }else{
    out<-newLG(cross = cross, markerList = chr.ord, keep.rf = T)
    return(out)
  }
}
