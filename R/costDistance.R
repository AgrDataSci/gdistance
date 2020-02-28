#' Cost distance
#' 
#' Calculate the least-cost distance between points.
#' 
#' @name costDistance
#' @rdname costDistance-methods
#' @aliases costDistance
#' @aliases costDistance,TransitionLayer,Coords,missing-method
#' @aliases costDistance,TransitionLayer,Coords,Coords-method
#' @keywords spatial
#' 
#' @param x object of class \code{TransitionLayer}
#' @param fromCoords first set of point locations
#'  (of \code{SpatialPoints}, matrix or numeric class)
#' @param toCoords optional, second set of point 
#'  locations (of \code{SpatialPoints}, matrix or numeric class)
#' @return distance matrix (S3 class dist or matrix)
#' @details 
#' Cost units between cells are defined as the reciprocal 
#'  of the values in the transition matrix.
#'  
#' The function uses Dijkstra's algorithm, as implemented in 
#'  the igraph package.
#'  
#' A projection correction is needed for accuracy in the case 
#'  of grid data for a longlat raster (see function geoCorrection).
#'  
#' @references 
#' E.W. Dijkstra. 1959. A note on two problems 
#'  in connexion with graphs. Numerische Mathematik 1, 269 - 271.
#' @author Jacob van Etten
#' @seealso 
#'  \code{\link{geoCorrection}},
#'  \code{\link{accCost}},
#' @examples
#' library("raster")
#' # create a new raster and set all its values to unity.
#' r <- raster(nrows=18, ncols=36) 
#' r <- setValues(r,runif(ncell(r),0,1))
#' 
#' # create a Transition object from the raster
#' tr <- transition(r,function(x) 1/mean(x),8)
#' 
#' # asymmetric
#' ncf <- function(x) max(x) - x[1] + x[2]
#' tr2 <- transition(r,ncf,8, symm=FALSE)
#' 
#' # create two sets of coordinates
#' sP1 <- cbind(c(65,5,-65),c(55,35,-35))
#' sP2 <- cbind(c(50,15,-40),c(80,20,-5))
#' 
#' # from and to identical
#' costDistance(tr,sP1)
#' costDistance(tr2,sP1)
#' 
#' # from and to different
#' costDistance(tr,sP1,sP2)
#' costDistance(tr2,sP1,sP2)
#' 
#' @exportMethod costDistance
setGeneric("costDistance", 
           function(x, fromCoords, toCoords) standardGeneric("costDistance"))

setMethod("costDistance", 
          signature(x = "TransitionLayer", 
                    fromCoords = "Coords", 
                    toCoords = "Coords"), 
          def = function(x, fromCoords, toCoords)
{
  return(.cd(x, fromCoords, toCoords))
}
)

.cd <- function(x, fromCoords, toCoords)
{
  fromCoords <- .coordsToMatrix(fromCoords)
  toCoords <- .coordsToMatrix(toCoords)
  
  fromCells <- raster::cellFromXY(x, fromCoords)
  toCells <- raster::cellFromXY(x, toCoords)
  
  if(!all(!is.na(fromCells))){
    warning("some coordinates not found and omitted")
    fromCells <- fromCells[!is.na(fromCells)]
  }
  
  if(!all(!is.na(toCells))){
    warning("some coordinates not found and omitted")
    toCells <- toCells[!is.na(toCells)]
  }
  
  costDist <- matrix(NA, nrow=length(fromCoords[,1]),ncol=length(toCoords[,1]))
  rownames(costDist) <- rownames(fromCoords)
  colnames(costDist) <- rownames(toCoords)
  y <- transitionMatrix(x)
  if(isSymmetric(y)) {m <- "undirected"} else{m <- "directed"}
  adjacencyGraph <- igraph::graph.adjacency(y, mode=m, weighted=TRUE)
  
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  
  uniqueFromCells <- unique(fromCells)
  uniqueToCells <- unique(toCells)		
  shortestPaths <- igraph::shortest.paths(adjacencyGraph, 
                                          v=uniqueFromCells,
                                          to=uniqueToCells,
                                          mode="out", 
                                          algorithm="dijkstra")
  
  index1 <- match(fromCells,uniqueFromCells)
  index2 <- match(toCells,uniqueToCells)
  costDist[] <- shortestPaths[index1,index2]
  return(costDist)
}

setMethod("costDistance", 
          signature(x = "TransitionLayer", 
                    fromCoords = "Coords", 
                    toCoords = "missing"), 
          def = function(x, fromCoords)
{
  return(.cd2(x, fromCoords))
}
)

.cd2 <- function(x, fromCoords)
{
  fromCoords <- .coordsToMatrix(fromCoords)
  fromCells <- raster::cellFromXY(x, fromCoords)
  
  if(!all(!is.na(fromCells))){
    warning("some coordinates not found and omitted")
    fromCells <- fromCells[!is.na(fromCells)]
  }
  
  costDist <- matrix(NA, 
                     nrow=length(fromCoords[,1]),
                     ncol=length(fromCoords[,1]))
  rownames(costDist) <- rownames(fromCoords)
  colnames(costDist) <- rownames(fromCoords)
  
  if(isSymmetric(transitionMatrix(x))) {
    m <- "undirected"
  }else{
      m <- "directed"
  }
  
  adjacencyGraph <- igraph::graph.adjacency(transitionMatrix(x),
                                            mode=m, weighted=TRUE)
  
  E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight
  
  uniqueFromCells <- unique(fromCells)
  
  shortestPaths <- igraph::shortest.paths(adjacencyGraph,
                                          v=uniqueFromCells, 
                                          to=uniqueFromCells,
                                          mode="out")
  
  index <- match(fromCells,uniqueFromCells)
  costDist[] <- shortestPaths[index,index]
  
  if(m=="undirected") {
    costDist <- as.dist(costDist)
  }
  
  return(costDist)
}

# TO DO check if coordinate systems are equal.
# TO DO check if bounding box of coordinates falls inside bb of transition

