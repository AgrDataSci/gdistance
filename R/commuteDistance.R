#' Commute-time distance
#'  
#' Calculates the resistance distance between points.
#' 
#' @name commuteDistance
#' @rdname commuteDistance-methods
#' @aliases commuteDistance,TransitionLayer,Coords-method
#' @keywords spatial
#' 
#' @param x object of class \code{TransitionLayer}
#' @param coords point locations coordinates 
#' (of SpatialPoints, matrix or numeric class)
#' @return distance matrix (S3 class dist or matrix)
#' @details 
#' This function calculates the expected random-walk commute 
#'  time between nodes in a graph. It is defined as the effective 
#'  distance (resistance distance) between the selected nodes multiplied 
#'  by the volume of the graph, which is the sum of the conductance 
#'  weights of all the edges in the graph (Chandra et al. 1997). The 
#'  result represents the average number of steps that is needed to 
#'  commute between the nodes during a random walk.
#' 
#' The function implements the algorithm given by Fouss et al. (2007).
#' 
#' Before calculating commute-time distances from a \code{TransitionLayer} 
#'  object, see if you need to apply the function 
#'  \code{\link{geoCorrection}}
#' 
#' @references 
#' Chandra, A.K., Raghavan, P., Ruzzo, W.L., Smolensy, R. & Tiwari, P. 1996. 
#'  The electrical resistance of a graph captures its commute and cover times. 
#'  Computational Complexity, 6(4), 312-340.
#' 
#' Fouss, F., Pirotte, A., Renders, J.-M. & Saerens, M. 2007. 
#'  Random-walk computation of similarities between nodes of a 
#'  graph with application to collaborative recommendation. IEEE 
#'  Transactions on Knowledge and Data Engineering, 19(3), 355-369.
#' 
#' McRae, B.H. 2006. Isolation by resistance. Evolution 60(8), 1551-1561.
#' 
#' \url{http://www.circuitscape.org}
#' 
#' @author Jacob van Etten
#' @seealso \code{\link{geoCorrection}}
#' @examples
#' library("raster")
#' # Create a new raster and set all its values to unity.
#' r <- raster(nrows=18, ncols=36) 
#' r <- setValues(r,rep(1,ncell(raster)))
#' 
#' #Create a Transition object from the raster
#' tr <- transition(r, function(x) 1/mean(x),4)
#' 
#' # Create two sets of coordinates
#' library("sp")
#' sP1 <- SpatialPoints(cbind(c(65,5,-65),c(55,35,-35)))
#' sP2 <- SpatialPoints(cbind(c(50,15,-40),c(80,20,-5)))
#' 
#' #Calculate the resistance distance between the points
#' commuteDistance(tr, sP1)
#' 
#' @exportMethod commuteDistance
setGeneric("commuteDistance", 
           function(x, coords) standardGeneric("commuteDistance"))

setMethod("commuteDistance", 
          signature(x = "TransitionLayer", coords = "Coords"), 
          def = function(x, coords) 
{
  return(.rD(x, coords))
}
)

setMethod("commuteDistance", 
          signature(x = "TransitionLayer", coords = "Coords"),
          def = function(x, coords) 
{
  return(.rD(x, coords))
}
)

.rD <- function(x, coords){
  if(class(transitionMatrix(x)) != "dsCMatrix"){
    stop("symmetric transition matrix required",
         "(dsCMatrix) in TransitionLayer object x")}
  coords <- .coordsToMatrix(coords)
  
  rd <- matrix(NA,nrow=length(coords[,1]),ncol=length(coords[,1]))
  rownames(rd) <- rownames(coords)
  colnames(rd) <- rownames(coords)
  allFromCells <- cellFromXY(x, coords)
  
  if(!all(!is.na(allFromCells))){
    warning("some coordinates not found and omitted")
    allFromCells <- allFromCells[!is.na(allFromCells)]
  }
  
  x <- .transitionSolidify(x)		
  fromCells <- allFromCells[allFromCells %in% transitionCells(x)]
  if (length(fromCells) < length(allFromCells)) 
  {
    warning(length(fromCells)," out of ",length(allFromCells),
            " locations were found in the fully connected transition",
            " matrix. NAs introduced.")
  }
  else{}
  fromCells <- unique(allFromCells)
  
  
  Lr <- .Laplacian(x)
  n <- max(Lr@Dim)
  Lr <- Lr[-n,-n]
  #This should avoid too big floating points as "Voltage differences",
  # but give a number that can still be divided by n
  C <- 1e-300 * n 
  Lplus <- matrix(ncol=length(fromCells),nrow=length(fromCells))
  index <- match(fromCells,transitionCells(x))
  #Lr <- Cholesky(Lr)
  for (i in 1:length(fromCells))
  {
    ei <- matrix(-C/n, ncol=1, nrow=n-1)
    ei[index[i],] <- C-(C/n)
    xi <- solve(Lr,ei)
    xi <- as.vector(xi)
    Lplusallrows <- c(xi-sum(xi/n),(sum(xi)/n)) 
    Lplus[,i] <- as.vector(Lplusallrows)[index]
  }
  Lplus <- Lplus / C
  rdSS <- (-2*Lplus + matrix(diag(Lplus),nrow=length(fromCells),
                             ncol=length(fromCells)) 
           + t(matrix(diag(Lplus),nrow=length(fromCells),
                      ncol=length(fromCells)))) 
  
  Volume <- sum(transitionMatrix(x))
  rdSS <- rdSS * Volume
  
  index1 <- which(allFromCells %in% fromCells)
  index2 <- match(allFromCells[allFromCells %in% fromCells],fromCells)
  rd[index1,index1] <- rdSS[index2,index2]
  rd <- as.dist(rd)
  attr(rd, "method") <- "commute"
  return(rd)
}

# TO DO check if coordinate systems are equal.
# TO DO check if bounding box of coordinates falls inside bb of transition
# TO DO coordinates in same cell: distance = 0

