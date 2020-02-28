#' Probability of passage
#' 
#' Calculates for each cell the number of passages of a random-walk 
#'  or randomised shortest paths with given origin(s) and destination(s). 
#'  Either the total or the net number of passages can be calculated.
#'  In the case of multiple origins or destinations, each receives equal 
#'  weight. 
#' 
#' @name passage
#' @aliases passage
#' @aliases passage,TransitionLayer,Coords,Coords,missing-method
#' @aliases passage,TransitionLayer,Coords,Coords,numeric-method
#' @aliases passage,TransitionLayer,RasterLayer,RasterLayer,missing-method
#' @aliases passage,TransitionLayer,RasterLayer,RasterLayer,numeric-method
#' @keywords spatial
#' @keywords methods
#' @param x Object of class \code{Transition*}
#' @param origin \code{SpatialPoints}, matrix or numeric 
#'  object with coordinates or RasterLayer object with origin 
#'  cells set to TRUE
#' @param goal \code{SpatialPoints}, matrix or numeric object 
#'  with coordinates or RasterLayer object with origin 
#'  cells set to TRUE
#' @param theta If zero or missing, a random walk results.
#'  If a numeric value 0 < theta < 20 is given, randomised 
#'  shortest paths are calculated; theta is the degree
#'  from which the path randomly deviates from the shortest path 
#' @param ... Additional arguments: 
#'   totalNet ("total" or "net"), and output ("RasterLayer" or "Transition")
#' @return RasterLayer or Transition object, depending on the output argument
#' @details 
#' The net number of passages between i and j is defined as: 
#'  abs( passages from i to j - passages from j to i ).
#'  
#' Defaults for additional argument \code{totalNet} is "net" 
#'  and for \code{output} it is "RasterLayer".
#' 
#' Random walk requires a symmetric transition matrix.
#' 
#' @references 
#' McRae B.H., B.G. Dickson, and T. Keitt. 2008. 
#'  Using circuit theory to model connectivity in ecology, 
#'  evolution, and conservation. 
#'  \emph{Ecology} 89:2712-2724.
#' 
#' Saerens M., L. Yen, F. Fouss, and Y. Achbany. 2009. 
#'  Randomized shortest-path problems: two related models. 
#'  \emph{Neural Computation}, 21(8):2363-2404.
#' 
#' @author Jacob van Etten. Implementation of randomised shortest 
#'  paths based on Matlab code by Marco Saerens
#'  
#' @seealso 
#' \code{\link[gdistance]{commuteDistance}}, 
#' \code{\link[gdistance]{pathInc}}
#' 
#' @examples 
#' #create a new raster and set all its values to unity.
#' raster <- raster(nrows=18, ncols=36) 
#' raster <- setValues(raster,rep(1,ncell(raster)))
#' 
#' #create a Transition object from the raster
#' tr <- transition(raster,mean,4)
#' trC <- geoCorrection(tr, type="c", scl=TRUE)
#' trR <- geoCorrection(tr, type="r", scl=TRUE)
#' 
#' #create two coordinates
#' sP1 <- SpatialPoints(cbind(-105,55))
#' sP2 <- SpatialPoints(cbind(105,-55))
#' 
#' #randomised shortest paths with theta = 2
#' rSPraster <- passage(trC, sP1, sP2, 2)
#' plot(rSPraster)
#' points(sP1)
#' points(sP2)
#' 
#' #randomised shortest paths with theta = 0.05
#' rSPraster <- passage(trC, sP1, sP2, 0.05)
#' plot(rSPraster) 
#' points(sP1)
#' points(sP2)
#' 
#' #randomised shortest paths with theta = 0.05
#' #and total
#' rSPraster <- passage(trC, sP1, sP2, 0.05, totalNet = "total")
#' plot(rSPraster) 
#' points(sP1)
#' points(sP2)
#' 
#' #random walk
#' rwraster <- passage(trR, sP1, sP2)
#' plot(rwraster)
#' points(sP1)
#' points(sP2)
#' 
#' @exportMethod passage
setGeneric("passage", function(x, origin, goal, theta, ...) {
  standardGeneric("passage")
}
)

setMethod("passage", 
          signature(x = "TransitionLayer", 
                    origin = "Coords", 
                    goal = "Coords", 
                    theta="missing"), 
          def = function(x, origin, goal, 
                         totalNet="net", 
                         output="RasterLayer")
	{
		.checkInputsPassage(totalNet, output)
		
		origin <- .coordsToMatrix(origin)
		goal <- .coordsToMatrix(goal)
		
		if(totalNet=="net" & output=="RasterLayer")
		{
			x <- .transitionSolidify(x)
			tc <- transitionCells(x)
			cellnri <- cellFromXY(x, origin)
			cellnrj <- cellFromXY(x, goal)
			ci <- match(cellnri,tc)
			cj <- match(cellnrj,tc)
			result <- .flowMap(x, ci, cj, tc)
		}
		else{
		  stop("no method available -- try a low value of theta instead")
		}
		
		return(result)
	}
)

setMethod("passage", signature(x = "TransitionLayer", 
                               origin = "RasterLayer", 
                               goal = "RasterLayer", 
                               theta="missing"), 
          def = function(x, origin, goal, 
                         totalNet="net", 
                         output="RasterLayer")
	{
		.checkInputsPassage(totalNet, output)
		
		if(totalNet=="net" & output=="RasterLayer")
		{
			x <- .transitionSolidify(x)
			tc <- transitionCells(x)
			ci <- which(getValues(origin))
			cj <- which(getValues(goal))
			result <- .flowMap(x, ci, cj, tc)
		}
		else{
		  stop("no method available -- try a low value of theta instead")
		}
            
		return(result)
	}
)

.checkInputsPassage <- function(totalNet, output)
{
	if(!(totalNet %in% c("total","net"))){
	  stop("totalNet should be either total or net")
	}
  
	if(!(output %in% c("RasterLayer","TransitionLayer"))){
	  stop("output should be either RasterLayer or TransitionLayer")
	}
  
} 

.flowMap <- function(x, indexOrigin, indexGoal, tc)
{
	L <- .Laplacian(x)
	Lr <- L[-dim(L)[1],-dim(L)[1]]
	A <- as(L,"lMatrix")
	A <- as(A,"dMatrix")
	n <- max(dim(Lr))
	Current <- .currentR(L, Lr, A, n, indexOrigin, indexGoal)
	result <- as(x,"RasterLayer")
	dataVector <- rep(0,times=ncell(result))
	dataVector[tc] <- Current
	result <- setValues(result, dataVector)
	return(result)
}

# Author: Jacob van Etten jacobvanetten@yahoo.com, 
# based on Matlab code by Marco Saerens
# IE University
# Date :  January 2010
# Version 1.0
# Licence GPL v3

setMethod("passage", signature(x = "TransitionLayer", 
                               origin = "Coords", 
                               goal = "Coords", 
                               theta="numeric"), 
          def = function(x, origin, goal, theta, 
                         totalNet="net",
                         output="RasterLayer")
	{
		cellnri <- cellFromXY(x, origin)
		cellnrj <- cellFromXY(x, goal)
		x <- .transitionSolidify(x)
		tc <- transitionCells(x)

		ci <- match(cellnri,tc)
		cj <- match(cellnrj,tc)
		
		result <- .randomShPaths(x, ci, cj, theta, 
		                         tc, totalNet, output)
		return(result)
	}
)

setMethod("passage", signature(x = "TransitionLayer", 
                               origin = "RasterLayer", 
                               goal = "RasterLayer", 
                               theta="numeric"), 
          def = function(x, origin, goal, theta, 
                         totalNet="net", output="RasterLayer")
	{
		#check if Transition and RasterLayers coincide
		ci <- which(getValues(origin))
		cj <- which(getValues(goal))
		x <- .transitionSolidify(x)
		tc <- transitionCells(x)
		result <- .randomShPaths(x, ci, cj, theta, tc, totalNet, output)
		return(result)
	}
)

.randomShPaths <- function(x, ci, cj, theta, tc, totalNet, output)
{
	if(theta < 0 | theta > 20 ) {
	  stop("theta value out of range (between 0 and 20)")
	}
  
	
	tr <- transitionMatrix(x,inflate=FALSE)
	
	trR <- tr
	trR@x <- 1 / trR@x 
	nr <- dim(tr)[1] 
	Id <- Diagonal(nr) 
	rs <- rowSums(tr)
	rs[rs>0] <- 1/rs[rs>0]
	P <- tr * rs

	W <- trR
	#zero values are not relevant because of next step exp(-theta * trR@x) 
	W@x <- exp(-theta * trR@x) 
	W <- W * P 

	return(.probPass(x, Id, W, nr, ci, cj, tc, totalNet, output))
}
	
.probPass <- function(x, Id, W, nr, ci, cj, tc, totalNet, output)
{
	nc <- ncell(x)
	
	Ij <- Diagonal(nr)
	Ij[cbind(cj,cj)] <- 1 - 1 / length(cj)
	Wj <- Ij %*% W
	
	ei <- rep(0,times=nr)
	ei[ci] <- 1 / length(ci)
	
	ej <- rep(0,times=nr)

	ej[cj] <- 1 / length(cj)
	
	IdMinusWj <- as((Id - Wj), "dgCMatrix")
	
	zci <- solve(t(IdMinusWj),ei)
	zcj <- solve(IdMinusWj, ej)
	zcij <- sum(ei*zcj)
	
	if(zcij < 1e-300)
	{
		if(output == "RasterLayer")	
		{
			result <- as(x,"RasterLayer")
			result[] <- rep(0,times=nc)
		}	
		if(output == "TransitionLayer")
		{
			result <- x
			transitionMatrix(result) <- Matrix(0, nc, nc)
			result@transitionCells <- 1:nc
		}
	}
	
	else
	{
		# Computation of the cost dij between node i and node j
		# dij <- (t(zci) %*% (trR * Wj) %*% zcj) / zcij
	
		# Computation of the matrix N, containing the number of passages through
		# each arc
		N <- (Diagonal(nr, as.vector(zci)) %*% Wj %*% Diagonal(nr, as.vector(zcj))) / zcij
	
		#N is here the NET number of passages, like McRae-random walk

		if(output == "RasterLayer")
		{
			if(totalNet == "total")
			{
				# Computation of the vector n, containing the number of visits in
				# each node
				n <- pmax(rowSums(N),colSums(N)) #not efficient but effective
			}
			# Computation of the matrix Pr, containing the transition
			# probabilities
			#rn <- rep(0, times=length(n))
			#rn[n>0] <- 1 / n[n>0]
			#Pr <- N * rn
			#Pr <- N * (1 / n)
	
			#net visits
			if(totalNet == "net")
			{
				nNet <- abs(skewpart(N))
				n <- pmax(rowSums(nNet),colSums(nNet))
				n[c(ci,cj)] <- 2 * n[c(ci,cj)]
			}
			result <- as(x,"RasterLayer")
			dataVector <- rep(NA,times=nc)
			dataVector[tc] <- n
			result <- setValues(result, dataVector)
		}
	
		if(output == "TransitionLayer")
		{
			result <- x
			if(totalNet == "total")
			{
				tr <- Matrix(0, nc, nc)
				tr[tc,tc] <- N
			}
			if(totalNet == "net")
			{
				nNet <- skewpart(N) * 2
				nNet@x[nNet@x<0] <- 0
				tr <- Matrix(0, nc, nc)
				tr[tc,tc] <- nNet
			}
			
			transitionMatrix(result) <- tr
			result@transitionCells <- 1:nc
			
		}
	}
	return(result)
}

