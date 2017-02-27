#' @title Order markers within linkage groups
#'
#' @description
#' \code{tspOrder} Apply a traveling salesperson problem solver to find the
#' shortest path through the recombination fraction matrix. Requires that
#' the TSP package is installed. Optimal performance is provided by the
#' "Concorde" method; however, this requires independent installation of
#' the Concorde program. See details.
#'
#' @param cross The QTL cross object.
#' @param method The solve_TSP method to employ. We highly encourage using method = "condcorde".
#' Simulations show that this method outperforms all others significantly. See details.
#' @param hamiltonian Logical, should a hamiltonian circuit be enforced?
#' @param concorde_path Required if method = "concorde". The directory containing
#' concorde executables.
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
#' The best performance seems to result from using hamiltonian = T and method = "concorde".
#'
#' We recommend using the concorde algorithm. To do this, install the concorde program:
#' http://www.math.uwaterloo.ca/tsp/concorde/downloads/codes/src/co031219.tgz.
#' If on a mac, this is not super easy.
#' See https://qmha.wordpress.com/2015/08/20/installing-concorde-on-mac-os-x/
#' for details on the best way to do this.
#'
#' @return Either a cross object with reordered markers, or a named list
#' of markers in a new order.
#'
#' @examples
#' library(qtlTools)
#' data(fake.f2)
#' cross<-fake.f2
#' \dontrun{
#' fake.f2<-est.rf(fake.f2)
#' cross<-fake.f2
#' #Perturb the marker order and chromosome names
#' markerlist<-lapply(chrnames(cross), function(x) sample(markernames(cross, chr =x)))
#' names(markerlist)<-as.character(chrnames(cross))
#' cross2<-newLG(cross, markerList = markerlist)
#' library(TSP)
#' plot.rf(cross2)
#' cross3<-cross3<-tspOrder(cross = cross2,
#'   hamiltonian = T,
#'   method="nn") # change to your path
#' cross3<-cross3<-tspOrder(cross = cross2,
#'   hamiltonian = T,
#'   method="concorde",
#'   concorde_path = "/Users/John/Documents/concorde/TSP") # change to your path
#' plot.rf(cross3)
#' }
#' @import qtl
#' @export
tspOrder<-function(cross,
                   method = "concorde",
                   hamiltonian = TRUE,
                   concorde_path = NULL,
                   return = "cross"){

  newLG<-function(cross, markerList){
    if(any(is.null(names(markerList)))){
      names(markerList)<-as.character(1:length(markerList))
    }
    # drop markers not in markerList
    newmars<-unlist(markerList)
    if(any(duplicated(newmars))){
      stop("duplicated markers found in markerList, all markers must be unique\n")
    }
    oldmars<-markernames(cross)
    cross <- drop.markers(cross, oldmars[!oldmars %in% newmars])



    # pull out rf matrix before reformatting cross
    n.mar <- nmar(cross)
    tot.mar <- totmar(cross)
    has.rf<-ifelse("rf" %in% names(cross), TRUE, FALSE)
    if(has.rf){
      rf <- cross$rf
      lod <- rf
      lod[lower.tri(rf)] <- t(rf)[lower.tri(rf)]
      rf[upper.tri(rf)] <- t(rf)[upper.tri(rf)]
      diagrf <- diag(rf)
      diag(lod) <- 0
      if(ncol(rf) != tot.mar)
        stop("dimension of recombination fractions inconsistent with no. markers in cross.")
      marnam <- colnames(rf)
    }

    chrstart <- rep(names(cross$geno), n.mar)

    # clean the cross
    cross <- clean(cross)
    crosstype <- class(cross)[1]
    g <- pull.geno(cross)


    cross$geno <- vector("list", length(markerList))
    names(cross$geno) <- names(markerList)

    for(i in names(markerList)) {
      cross$geno[[i]]$data <- g[,markerList[[i]],drop=FALSE]

      cross$geno[[i]]$map <- seq(0, by=10, length=length(markerList[[i]]))
      if(crosstype=="4way") {
        cross$geno[[i]]$map <- rbind(cross$geno[[i]]$map,
                                     cross$geno[[i]]$map)
        colnames(cross$geno[[i]]$map) <- markerList[[i]]
      }else{
        names(cross$geno[[i]]$map) <- markerList[[i]]
      }
      class(cross$geno[[i]]) <- "A"
    }

    if(has.rf){
      mname <- markernames(cross)
      m <- match(mname, marnam)
      rf <- rf[m,m]
      lod <- lod[m,m]
      rf[upper.tri(rf)] <- lod[upper.tri(lod)]
      diag(rf) <- diagrf[m]
      cross$rf <- rf
    }
    if("X" %in% chrnames(cross)){
      class(cross$geno[["X"]]) <- "X"
    }
    return(cross)
  }

  if(!requireNamespace("TSP", quietly = TRUE)){
    stop("install the TSP package to use tspOrder\n")
  }else{
    requireNamespace("TSP", quietly = TRUE)
  }
  if(method == "concorde"){
    if(is.null(concorde_path)){
      stop("if method = concorde, concorde_path must specify directory of program\n")
    }else{
      concorde_path(concorde_path)
    }
  }
  rf<-data.matrix(pull.rf(cross, what = "rf"))
  class(rf)<-"matrix"
  diag(rf)<-0

  markerList<-lapply(chrnames(cross), function(x) markernames(cross, chr = x))
  names(markerList)<-chrnames(cross)

  chr.ord<-lapply(names(markerList),function(x){
    mnames<-markerList[[x]]
    totsp<-rf[mnames,mnames]
    chr.tsp<-TSP(totsp)
    if(hamiltonian){
      tsp <- insert_dummy(chr.tsp, label = "cut")
      tour <- solve_TSP(tsp, method=method)
      path <- cut_tour(tour, "cut")
    }else{
      path <- solve_TSP(chr.tsp, method=method)
    }
    ord<-labels(path)
    return(ord)
  })
  names(chr.ord)<-names(markerList)
  if(return != "cross"){
    return(chr.ord)
  }else{
    out<-newLG(cross = cross, markerList = chr.ord)
    return(out)
  }
}
