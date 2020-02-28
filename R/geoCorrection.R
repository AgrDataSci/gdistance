#' Geographic Correction
#' 
#' Correct 'TransitionLayer' objects taking into account 
#'  local distances
#'  
#' @name geoCorrection
#' @aliases geoCorrection
#' @aliases geoCorrection,TransitionLayer,character-method
#' @aliases geoCorrection,TransitionLayer,missing-method
#' @keywords spatial
#' @keywords methods
#' 
#' @param x object of class 'Transition*'
#' @param type type of correction: "c", "r", or missing 
#'  (only required for lonlat, see Details)
#' @param ... additional arguments passed to methods. \code{multpl} with 
#'  correction factor (\code{TRUE}) or corrected values (\code{FALSE}, the default).
#'  \code{scl} scale the correction values (default is \code{FALSE})
#' @return a 'Transition*' object
#' @details 
#' 
#' Geographic correction is necessary for all objects of the 
#'  class Transition that are either: 
#'  (1) based on a grid in a geographic (lonlat) projection 
#'   and covering a large area; 
#'  (2) made with directions > 4.
#' 
#' The function will correct for map distortion, as well as for diagonal 
#'  connections between grid cells (which cover a longer distance than vertical 
#'  or horizontal connections).
#' 
#' When working with lonlat grids, users should also anticipate whether they 
#'  will use methods based on either least-cost or random walks, and set the 
#'  type argument accordingly. In the case of least-cost distances, the 
#'  correction is only done in East-West direction. In the case of random walks 
#'  there should be an additional correction which reduces the conductance in
#'  North-South direction (\code{type="r"}).
#' 
#' The correction is done by dividing conductance values by the inter-cell 
#'  distance. Distances are calculated as great-circle distances for lonlat 
#'  grids (see function \code{[raster]{isLonLat}}) and Euclidean distances for all other grids.
#' 
#' In the case of randomised shortest paths, the need for correction is 
#' somewhat in between these two correction methods. We have not developed 
#' an analytical solution for this problem yet. With very low values for theta,
#' you may choose the correction for random walks, and for high values the one for 
#' least-cost paths. Users who want to work with intermediate values of theta are 
#' encouraged to experiment with different solutions.
#' 
#' The values are scaled to get values near 1 if the argument \code{scl} is set to 
#'  TRUE. This is desirable for subsequent calculations involving random walk
#'  calculations. Values are divided by the W-E inter-cell distance (at the 
#'  centre of the grid).
#' 
#' @author Jacob van Etten
#' 
#' @examples
#' library("raster")
#' r <- raster(ncol=36,nrow=18)
#' r <- setValues(r,rep(1,times=ncell(r)))
#' tr <- transition(r, mean, directions=8)
#' 
#' # directly
#' tr1 <- geoCorrection(tr, type="c", multpl=FALSE)
#' tr1
#' 
#' # the same, but with a separate correction matrix
#' trCorr <- geoCorrection(tr, type="c", multpl=TRUE)
#' tr2 <- tr * trCorr
#' tr2
#' @exportMethod geoCorrection
setGeneric("geoCorrection", function(x, type, ...) {
  standardGeneric("geoCorrection")
  }
)


setMethod("geoCorrection", 
          signature(x = "TransitionLayer", type="missing"), 
          def = function(x, multpl=FALSE, scl=FALSE)
{
  return(geoCorrection(x, type="c", multpl, scl))
}
)


setMethod("geoCorrection", 
          signature(x = "TransitionLayer", type="character"),
          def = function(x, type, multpl=FALSE, scl=FALSE)
{
  scaleValue <- 1
  if(scl)
  {
    midpoint <- c(mean(c(xmin(x),xmax(x))),mean(c(ymin(x),ymax(x))))
    scaleValue <- pointDistance(midpoint,midpoint+c(xres(x),0), longlat=isLonLat(x)) 
  }
  if(isLonLat(x))
  {
    if (type != "c" & type != "r"){
      stop("type can only be c or r")
    }
    
    # if (type == "r" & matrixValues(x) != "conductance"){
    #   stop("matrix of Transition object must have conductance values")
    # }
    
    adj <- adjacencyFromTransition(x)
    correction <- cbind(raster::xyFromCell(x,adj[,1]),
                        raster::xyFromCell(x,adj[,2]))
    if(matrixValues(x) == "conductance") {
      correctionValues <- 1/(raster::pointDistance(correction[,1:2],
                                                   correction[,3:4],
                                                   longlat=TRUE)/scaleValue)
    }
    
    if(matrixValues(x) == "resistance") {
      correctionValues <- raster::pointDistance(correction[,1:2],
                                                correction[,3:4],
                                                longlat=TRUE)/scaleValue
    }
    
    if (type=="r")
    {
      rows <- raster::rowFromCell(x,adj[,1]) != raster::rowFromCell(x,adj[,2])
      #low near the poles
      if(matrixValues(x) == "conductance") {
        corrFactor <- cos((pi/180) * rowMeans(cbind(correction[rows,2],
                                                    correction[rows,4])))
      } 
      
      #high near the poles
      if(matrixValues(x) == "resistance") {
        corrFactor <- 1 / (cos((pi/180) * rowMeans(cbind(correction[rows,2],
                                                         correction[rows,4]))))
      } 
      
      # makes conductance lower in N-S direction towards the poles
      correctionValues[rows] <- correctionValues[rows] * corrFactor
    }
  } else {
    adj <- adjacencyFromTransition(x)
    correction <- cbind(raster::xyFromCell(x,adj[,1]),
                        raster::xyFromCell(x,adj[,2]))
    if(matrixValues(x) == "conductance") {
      correctionValues <- 1/(raster::pointDistance(correction[,1:2],
                                                   correction[,3:4],
                                                   longlat=FALSE)/scaleValue)
    }
    
    if(matrixValues(x) == "resistance") {
      correctionValues <- raster::pointDistance(correction[,1:2],
                                                correction[,3:4],
                                                longlat=FALSE)/scaleValue
    }
    
  }
  i <- as.integer(adj[,1] - 1)
  j <- as.integer(adj[,2] - 1)
  xv <- as.vector(correctionValues) #check for Inf values!
  dims <- ncell(x)
  correctionMatrix <- new("dgTMatrix", i = i, j = j, x = xv, 
                          Dim = as.integer(c(dims,dims)))
  correctionMatrix <- (methods::as(correctionMatrix,"sparseMatrix"))
  if(class(transitionMatrix(x)) == "dsCMatrix"){
    correctionMatrix <- forceSymmetric(correctionMatrix)
  }  #isSymmetric?
  if(!multpl) 
  {
    transitionCorrected <- correctionMatrix * transitionMatrix(x)
    transitionMatrix(x) <- transitionCorrected
    return(x)
  }	
  if(multpl)
  {
    transitionMatrix(x) <- correctionMatrix
    return(x)
  }
}
)

