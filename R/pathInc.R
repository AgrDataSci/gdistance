#' Incidence of paths from a common origin: overlap and non-overlap
#'
#' Calculate the overlap and non-overlap of paths departing from a common origin.
#'  Two algorithms are available: random walk and randomised shortest paths.
#'
#' @name pathInc
#' @aliases pathInc
#' @aliases pathInc,TransitionLayer,Coords,Coords,missing,missing,missing-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,missing,numeric,missing-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,missing,missing,Transition-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,missing,numeric,Transition-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,Coords,missing,missing-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,Coords,numeric,missing-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,Coords,missing,Transition-method
#' @aliases pathInc,TransitionLayer,Coords,Coords,Coords,numeric,Transition-method
#' @keywords spatial
#'
#' @param x transition matrix of class \code{Transition*}
#' @param origin coordinates of the origin (one point location,
#'  SpatialPoints, matrix or numeric class)
#' @param from coordinates of the destinations
#'  (SpatialPoints, matrix or numeric class)
#' @param to second set of coordinates of the destinations (can be missing)
#' @param theta value > 0 and < 20 (randomised shortest paths) or
#'  missing (random walk)
#' @param weight matrix – Reciprocals of the non-zero values are used as weights.
#' If missing, reciprocals of the transition matrix are used.
#' @param ... additional arguments passed to methods. See Details.
#' @return list of dist objects or list of matrices
#' @details
#' This is a complex wrapper function that calculates to what extent dispersal
#' routes overlap or do not overlap.
#'
#' First, the function calculates the trajectories for all "from" and "to"
#' locations, starting from a single "origin" location. These trajectories can
#' either be based on random walks or randomised shortest paths (giving a value
#' to theta).
#'
#' Then, for all unique pairs of trajectories, it calculates the extent to
#' which these trajectories overlap or diverge. These values are given back
#' to the user as a list of (distance) matrices.
#'
#' If only "from" coordinates are given, the function calculates symmetric
#' distance matrices for all combinations of points in "from". If both "from"
#' and "to" coordinates are given, the function calculates a matrix of values
#' with rows for all locations in "from" and columns for all locations in "to".
#'
#' Overlap is currently calculated as the minimum values of each pair of
#' trajectories compared. Non-overlap uses the following formula:
#' Nonoverlap = max(0,max(a,b)*(1-min(a,b))-min(a,b))
#' (see van Etten and Hijmans 2010). See the last example to learn
#' how to use an alternative function.
#'
#' @references
#' McRae B.H., B.G. Dickson, and T. Keitt. 2008.
#'  Using circuit theory to model connectivity in ecology,
#'  evolution, and conservation. Ecology 89:2712-2724.
#'
#' Saerens M., L. Yen, F. Fouss, and Y. Achbany. 2009.
#'  Randomized shortest-path problems: two related models.
#'  Neural Computation, 21(8):2363-2404.
#'
#' van Etten, J., and R.J. Hijmans. 2010. A geospatial modelling
#'  approach integrating archaeobotany and genetics to trace the
#'  origin and dispersal of domesticated plants. PLoS ONE 5(8): e12060.
#'
#' @author Jacob van Etten. Implementation of randomised shortest paths
#'  based on Matlab code by Marco Saerens.
#' @examples
#' library("raster")
#' library("sp")
#'
#' # Create TransitionLayer
#' r <- raster(ncol=36,nrow=18)
#' r <- setValues(r,rep(1,times=ncell(r)))
#' tr <- transition(r,mean,directions=4)
#'
#' # Two different types of correction are required
#' trR <- geoCorrection(tr, type="r", multpl=FALSE)
#' trC <- geoCorrection(tr, type="c", multpl=FALSE)
#'
#' # Create TransitionStack
#' ts <- stack(trR, trR)
#'
#' # Points for origin and coordinates between which to calculate path (non)overlaps
#' sP0 <- SpatialPoints(cbind(0,0))
#' sP1 <- SpatialPoints(cbind(c(65,5,-65),c(-55,35,-35)))
#'
#' # Randomised shortest paths
#' # rescaling is needed: exp(-theta * trC) should give reasonable values
#' # divide by median of the non-zero values
#' trC <- trC / median(transitionMatrix(trC)@x)
#' pathInc(trC, origin=sP0, from=sP1, theta=2)
#'
#' # Random walk
#' pathInc(trR, origin=sP0, from=sP1)
#'
#' # TransitionStack as weights
#' pathInc(trR, origin=sP0, from=sP1, weight=ts)
#'
#' # Demonstrate use of an alternative function
#' #The current default is to take the minimum of each pair of layers
#'
#' altoverlap <- function(a, b)
#' {
#'   aV <- as.vector(a[,rep(1:ncol(a), each=ncol(b))])
#'   bV <- as.vector(b[,rep(1:ncol(b), times=ncol(a))])
#'   result <- matrix(aV * bV, nrow = nrow(a), ncol=ncol(a)*ncol(b))
#'   return(result)
#' }
#'
#' pathInc(trR, origin=sP0, from=sP1, weight=ts, functions=list(altoverlap))
#'
#' @exportMethod pathInc
setGeneric("pathInc", function(x, origin, from, to, theta, weight, ...) {
  standardGeneric("pathInc")
  })

# to = "missing"

setMethod("pathInc", signature(x = "TransitionLayer",
                               origin = "Coords",
                               from = "Coords",
                               to = "missing",
                               theta="missing",
                               weight="missing"),
          def = function(x, origin, from,
                         functions=list(overlap,nonoverlap))
	{

		preparedMatrix <- .prepareMatrix(x, 0)
		preparedIndex <- .prepareIndex(preparedMatrix$transition, origin, from)
		prepared <- c(preparedMatrix,preparedIndex)
		Intermediate <- .randomWalk(prepared)
		result <- .finishFlow(prepared, Intermediate, functions)
		return(result)
	}
)

setMethod("pathInc", signature(x = "TransitionLayer",
                               origin = "Coords",
                               from = "Coords",
                               to = "missing", theta="numeric",
                               weight="missing"),
          def = function(x, origin, from, theta,
                         functions=list(overlap,nonoverlap))
	{
		if(theta < 0 | theta > 20 ) {
		  stop("theta value out of range (between 0 and 20)")
		}

		preparedMatrix <- .prepareMatrix(x, 0)
		preparedIndex <- .prepareIndex(preparedMatrix$transition, origin, from)
		prepared <- c(preparedMatrix,preparedIndex)
		Intermediate <- .randomSP(prepared, theta)
		result <- .finishFlow(prepared, Intermediate, functions)
		return(result)
	}
)

setMethod("pathInc", signature(x = "TransitionLayer",
                               origin = "Coords", from = "Coords",
                               to = "missing", theta="missing",
                               weight="Transition"),
          def = function(x, origin, from, weight,
                         functions=list(overlap,nonoverlap))
	{
		preparedMatrix <- .prepareMatrix(x, weight)
		preparedIndex <- .prepareIndex(preparedMatrix$transition, origin, from)
		prepared <- c(preparedMatrix,preparedIndex)
		Intermediate <- .randomWalk(prepared)
		result <- .finishFlow(prepared, Intermediate, functions)
		return(result)
	}
)

setMethod("pathInc", signature(x = "TransitionLayer",
                               origin = "Coords",
                               from = "Coords",
                               to = "missing", theta="numeric",
                               weight="Transition"),
          def = function(x, origin, from, theta, weight,
                         functions=list(overlap,nonoverlap))
	{
		if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
		preparedMatrix <- .prepareMatrix(x, weight)
		preparedIndex <- .prepareIndex(preparedMatrix$transition, origin, from)
		prepared <- c(preparedMatrix,preparedIndex)
		Intermediate <- .randomSP(prepared, theta)
		result <- .finishFlow(prepared, Intermediate, functions)
		return(result)
	}
)

# to = "Coords"

setMethod("pathInc", signature(x = "TransitionLayer",
                               origin = "Coords", from = "Coords",
                               to = "Coords", theta="missing",
                               weight="missing"),
          def = function(x, origin, from, to,
                         functions=list(overlap,nonoverlap))
	{
		preparedMatrix <- .prepareMatrix(x, 0)
		preparedIndexFrom <- .prepareIndex(preparedMatrix$transition, origin, from)
		preparedIndexTo <- .prepareIndex(x, origin, to)
		IntermediateFrom <- .randomWalk(c(preparedMatrix,preparedIndexFrom))
		IntermediateTo <- .randomWalk(c(preparedMatrix,preparedIndexTo))
		result <- .finishFlowFromTo(preparedIndexFrom,
		                            preparedIndexTo,
		                            IntermediateFrom,
		                            IntermediateTo,
		                            functions)
		return(result)
	}
)

setMethod("pathInc", signature(x = "TransitionLayer", origin = "Coords", from = "Coords",
	to = "Coords", theta="numeric", weight="missing"),
	def = function(x, origin, from, to, theta, functions=list(overlap,nonoverlap))
	{
		if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
		preparedMatrix <- .prepareMatrix(x, 0)
		preparedIndexFrom <- .prepareIndex(preparedMatrix$transition, origin, from)
		preparedIndexTo <- .prepareIndex(preparedMatrix$transition, origin, to)
		preparedFrom <- c(preparedMatrix, preparedIndexFrom)
		preparedTo <- c(preparedMatrix, preparedIndexTo)
		IntermediateFrom <- .randomSP(preparedFrom, theta)
		IntermediateTo <- .randomSP(preparedTo, theta)
		result <- .finishFlowFromTo(preparedIndexFrom,
		                            preparedIndexTo,
		                            IntermediateFrom,
		                            IntermediateTo,
		                            functions)
		return(result)
	}
)

setMethod("pathInc", signature(x = "TransitionLayer", origin = "Coords", from = "Coords",
	to = "Coords", theta="missing", weight="Transition"),
	def = function(x, origin, from, to, weight, functions=list(overlap,nonoverlap))
	{
		preparedMatrix <- .prepareMatrix(x, weight)
		preparedIndexFrom <- .prepareIndex(preparedMatrix$transition, origin, from)
		preparedIndexTo <- .prepareIndex(preparedMatrix$transition, origin, to)
		preparedFrom <- c(preparedMatrix, preparedIndexFrom)
		preparedTo <- c(preparedMatrix, preparedIndexTo)
		IntermediateFrom <- .randomWalk(preparedFrom)
		IntermediateTo <- .randomWalk(preparedTo)
		result <- .finishFlowFromTo(preparedIndexFrom,
		                            preparedIndexTo,
		                            IntermediateFrom,
		                            IntermediateTo,
		                            functions)
		return(result)
	}
)

setMethod("pathInc", signature(x = "TransitionLayer", origin = "Coords", from = "Coords",
	to = "Coords", theta="numeric", weight="Transition"),
	def = function(x, origin, from, to, theta, weight, functions=list(overlap,nonoverlap))
	{
		if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
		preparedMatrix <- .prepareMatrix(x, weight)
		preparedIndexFrom <- .prepareIndex(preparedMatrix$transition, origin, from)
		preparedIndexTo <- .prepareIndex(preparedMatrix$transition, origin, to)
		preparedFrom <- c(preparedMatrix, preparedIndexFrom)
		preparedTo <- c(preparedMatrix, preparedIndexTo)
		IntermediateFrom <- .randomSP(preparedFrom, theta)
		IntermediateTo <- .randomSP(preparedTo, theta)
		result <- .finishFlowFromTo(preparedIndexFrom,
		                            preparedIndexTo,
		                            IntermediateFrom,
		                            IntermediateTo,
		                            functions)
		return(result)
	}
)


#this function prepares the transition matrix
#it removes non-zero rows/columns
#it also prepares the weight matrix converts its matrix values to resistance if needed
.prepareMatrix <- function(x, weight)
{
	x <- .transitionSolidify(x)

	A <- as(transitionMatrix(x,inflate=FALSE),"lMatrix")
	A <- as(A,"dMatrix")
	AIndex <- as(as(A, "generalMatrix"), "TsparseMatrix")
	index1 <- cbind(transitionCells(x)[as.integer(AIndex@i+1)],
	                transitionCells(x)[as.integer(AIndex@j+1)])
	#TODO use adjacencyFromTransition()
	index2 <- cbind(as.integer(AIndex@i+1),
	                as.integer(AIndex@j+1))
	#if symmetric? index <- index[index[,1] < index[,2],]

	Size <- nrow(index1)

	if(is(weight, "numeric")) {
		R <- 1/x[index2]
		R[R == Inf] <- 0
		R <- matrix(R, nrow=1)
	}

	if(is(weight, "TransitionLayer")) {
		if(matrixValues(weight) == "conductance") {
		  R <- 1/weight[index1]} else{R <- weight[index1]
		  }

		R[R == Inf] <- 0
		R <- matrix(R, nrow=1)
	}
	if(is(weight, "TransitionStack")) {
		R <- matrix(nrow = nlayers(weight), ncol = length(index1[, 1]))
		for(i in 1:nlayers(weight))
		{
			if(matrixValues(weight[[i]]) == "conductance"){
			  R[i,] <- 1/weight[[i]][index1]} else{R[i,] <- weight[[i]][index1]
			  }
		  }
		R[R == Inf] <- 0
	}

	result <- list(transition=x,
					index=index2,
					A=A,
					R=R,
					Size=Size)
	return(result)
}

# this function maps the coordinates to the rows/columns of the transition matrix
# the mapping is slightly complicated by the fact that not all cells
# in the original raster have rows/colums in the transition matrix
.prepareIndex <- function(x, origin, fromCoords)
{
	origin <- .coordsToMatrix(origin)
	fromCoords <- .coordsToMatrix(fromCoords)

	originCell <- cellFromXY(x, origin)

	if (!(originCell %in% transitionCells(x))) {
	  stop("the origin refers to a zero row/column in the transition matrix (unconnected)")
	}

	allFromCells <- cellFromXY(x, fromCoords)
	fromCells <- allFromCells[allFromCells %in% transitionCells(x)]
	if (length(fromCells) < length(allFromCells))
	{
		warning(length(fromCells)," out of ",
		        length(allFromCells[,1]),
		        " locations were found inside the transition matrix. NAs introduced.")
	}
	fromCells <- unique(fromCells)

	indexCoords <- match(fromCells,transitionCells(x))
	indexOrigin <- match(originCell,transitionCells(x))


	#TODO rename fromCells to cells simply, because function is used for "from" and "to"
	result <- list(transition=x,
			fromCoords=fromCoords,
			allFromCells=allFromCells,
			fromCells=fromCells,
			indexCoords=indexCoords,
			indexOrigin=indexOrigin)
	return(result)
}

.randomWalk <- function(prepared)
{
	x <- prepared$transition
	indexCoords <- prepared$indexCoords
	indexOrigin <- prepared$indexOrigin
	fromCells <- prepared$fromCells
	index <- prepared$index
	Size <- prepared$Size
	A <- prepared$A


	L <- .Laplacian(x)
	Lr <- L[-dim(L)[1],-dim(L)[1]]
	n <- max(Lr@Dim)
	Lr <- Cholesky(Lr)

	if(!(canProcessInMemory(x, length(fromCells)*10)))
	#depending on memory availability, currents are calculated in a
	#piecemeal fashion or all at once
	{
		filenm=rasterTmpFile()
		Flow <- raster(nrows=length(fromCells), ncols=Size)
		Flow <- writeStart(Flow, filenm, overwrite=TRUE)
		for(i in 1:(length(fromCells)))
		{
			matrixRow <- .currentM(L, Lr, A, n, indexOrigin, indexCoords[i], index)
			Flow <- writeValues(Flow, matrixRow, start=1)
		}
		Flow <- writeStop(Flow)
	}
	else
	{
		Flow <- matrix(nrow=length(fromCells), ncol=Size)
		for(i in 1:(length(fromCells)))
		{
			Flow[i,] <- .currentM(L, Lr, A, n, indexOrigin, indexCoords[i], index)
		}

		Flow <- raster(Flow)

	}
	return(Flow)
}


.randomSP <- function(prepared, theta)
{
	x <- prepared$transition
	cj <- prepared$indexCoords
	ci <- prepared$indexOrigin
	index <- prepared$index
	Size <- prepared$Size
	fromCells <- prepared$fromCells

	tr <- transitionMatrix(x, inflate=FALSE)
	tc <- transitionCells(x)

	trR <- tr
	trR@x <- 1 / tr@x

	nr <- dim(tr)[1]
	Id <- Diagonal(nr)
	rs <- rowSums(tr)
	rs[rs>0] <- 1/rs[rs>0]
	P <- tr * rs

	W <- trR
	#zero values are not relevant because of next step exp(-theta * trR@x)
	W@x <- exp(-theta * trR@x)
	W <- W * P

	if(!(canProcessInMemory(x, length(fromCells)*10)))
	#this does not take into account the exact memory needed for matrix solving...
	{
		filenm=rasterTmpFile()
		Flow <- raster(nrows=length(cj), ncols=Size)
		Flow <- writeStart(Flow, filenm, overwrite=TRUE)
		for(i in 1:(length(cj)))
		{
			matrixRow <- transitionMatrix(.probPass(x, Id, W, nr, ci, cj[i], tc,
			                                        totalNet="net",
			                                        output="TransitionLayer"),
			                              inflate=FALSE)[index]
			Flow <- writeValues(Flow, matrixRow, 1)
		}
		Flow <- writeStop(Flow)
	}
	else
	{
		Flow <- matrix(nrow=length(cj),ncol=Size)
		for(i in 1:(length(cj)))
		{
			Flow[i,] <- transitionMatrix(.probPass(x, Id, W, nr, ci, cj[i], tc,
			                                       totalNet="net",
			                                       output="TransitionLayer"),
			                             inflate=FALSE)[index]
		}

		Flow <- raster(Flow)

	}
	return(Flow)
}

.finishFlow <- function(prepared, Flow1, functions)
{

	fromCells <- prepared$fromCells
	allFromCells <- prepared$allFromCells
	fromCoords <- prepared$fromCoords

	R <- prepared$R
	if(is.matrix(R)){nR <- nrow(R)} else{nR <- 1}
	nF <- length(functions)
	n <- nR * nF

	result <- vector("list", n)
	f <- paste("function", 1:length(functions), sep="")
	l <- paste("layer", 1:nR, sep="")
	names(result) <- paste(rep(f, each=nR), rep(l, times=nF), sep="")
	resulti <- matrix(nrow=length(fromCells), ncol=length(fromCells))
	for(i in 1:n){result[[i]] <- resulti}

	tr1 <- blockSize(Flow1, n=1)
	for(i in 1:tr1$n)
	{
		chunk1 <- getValues(Flow1, row=tr1$row[i], nrows=tr1$nrows[i])
		#rows are cell transitions, columns are locations
		chunk1 <- matrix(chunk1,nrow=ncol(Flow1))
		for(j in i:tr1$n) #the crucial difference with the asymmetric case
		{
			chunk2 <- getValues(Flow1, row=tr1$row[j], nrows=tr1$nrows[j])
			#rows are cell transitions, columns are locations
			chunk2 <- matrix(chunk2,nrow=ncol(Flow1))
			index2 <- cbind(rep(tr1$row[i]:(tr1$row[i]+tr1$nrows[i]-1),each=tr1$nrows[j]),
				rep(tr1$row[j]:(tr1$row[j]+tr1$nrows[j]-1),times=tr1$nrows[i]))

			for(f in 1:nF)
			{
				Product <- functions[[f]](chunk1,chunk2)

				for(r in 1:nR)
				{
					index1 <- (f-1) * nR + r
					result[[index1]][index2] <- colSums(Product * R[r,])
				}
			}
		}
	}

	index1 <- which(allFromCells %in% fromCells)
	index2 <- match(allFromCells[allFromCells %in% fromCells], fromCells)

	for(i in 1:length(result))
	{
		ri <- matrix(nrow=length(allFromCells),ncol=length(allFromCells))
		rownames(ri) <- rownames(fromCoords)
		colnames(ri) <- rownames(fromCoords)
		resulti <- as.matrix(as.dist(t(result[[i]])))
		ri[index1,index1] <- resulti[index2,index2]
		result[[i]] <- as.dist(ri)
	}

	return(result)

}

.finishFlowFromTo <- function(preparedFrom, preparedTo, Flow1, Flow2, functions)
{

	fromCells <- preparedFrom$fromCells
	allFromCells <- preparedFrom$allFromCells
	fromCoords <- preparedFrom$fromCoords

	toCells <- preparedTo$fromCells
	allToCells <- preparedTo$allFromCells
	toCoords <- preparedTo$fromCoords

	R <- preparedFrom$R
	if(is.matrix(R)){nR <- nrow(R)} else{nR <- 1}
	nF <- length(functions)
	n <- nR * nF

	result <- vector("list", n)
	f <- paste("function", 1:length(functions), sep="")
	l <- paste("layer", 1:nR, sep="")
	names(result) <- paste(rep(f, each=nR), rep(l, times=nF), sep="")
	resulti <- matrix(nrow=length(fromCells), ncol=length(toCells))
	for(i in 1:n){result[[i]] <- resulti}

	tr1 <- blockSize(Flow1, n=1)
	tr2 <- blockSize(Flow2, n=1)
	for(i in 1:tr1$n)
	{
		chunk1 <- getValues(Flow1, row=tr1$row[i], nrows=tr1$nrows[i])
		#rows are cell transitions, columns are locations
		chunk1 <- matrix(chunk1,nrow=ncol(Flow1))
		for(j in 1:tr2$n)
		{
			chunk2 <- getValues(Flow2, row=tr2$row[j], nrows=tr2$nrows[j])
			#rows are cell transitions, columns are locations
			chunk2 <- matrix(chunk2,nrow=ncol(Flow2))
			index2 <- cbind(rep(tr1$row[i]:(tr1$row[i]+tr1$nrows[i]-1),each=tr2$nrows[j]),
				rep(tr2$row[j]:(tr2$row[j]+tr2$nrows[j]-1),times=tr1$nrows[i]))

			for(f in 1:nF)
			{
				Product <- functions[[f]](chunk1,chunk2)

				for(r in 1:nR)
				{
					index1 <- (f-1) * nR + r
					result[[index1]][index2] <- colSums(Product * R[r,])
				}
			}
		}
	}

	index1from <- which(allFromCells %in% fromCells)
	index2from <- match(allFromCells[allFromCells %in% fromCells], fromCells)
	index1to <- which(allToCells %in% toCells)
	index2to <- match(allToCells[allToCells %in% toCells], toCells)

	for(i in 1:length(result))
	{
		ri <- matrix(nrow=length(allFromCells),ncol=length(allToCells))
		rownames(ri) <- rownames(fromCoords)
		colnames(ri) <- rownames(toCoords)
#		resulti <- t(result[[i]])
		ri[index1from,index1to] <- resulti[index2from,index2to]
		result[[i]] <- as.dist(ri)
	}


	return(result)

}

#' Overlap and nonoverlap of trajectories
#'
#' Special functions to calculate the degree of overlap and
#'  nonoverlap between trajectories
#'
#' @aliases overlap
#' @aliases nonoverlap
#' @keywords spatial
#' @param a \code{matrix} object
#' @param b \code{matrix} object
#' @return A matrix object
#' @details
#' These functions are used by the \code{pathInc()} as defaults.
#' @note
#' The functions are provided here to assist the user in defining
#'  alternative measures of overlap / nonoverlap.
#' @author Jacob van Etten
#' @export
overlap <- function(a, b){
	aV <- as.vector(a[,rep(1:ncol(a), each=ncol(b))])
	bV <- as.vector(b[,rep(1:ncol(b), times=ncol(a))])
	result <- matrix(pmin(aV, bV), nrow = nrow(a), ncol=ncol(a)*ncol(b))
	return(result)
}


#' @export
nonoverlap <- function(a, b){
	aV <- as.vector(a[,rep(1:ncol(a), each=ncol(b))])
	bV <- as.vector(b[,rep(1:ncol(b), times=ncol(a))])
	result <- matrix(pmax(pmax(aV, bV) * (1-pmin(aV, bV)) - pmin(aV, bV), 0),
	                 nrow = nrow(a),
	                 ncol=ncol(a)*ncol(b))
	return(result)
}

#TO DO check vignette example
#TO DO from+to case in .Rd file
#TO DO check all cases from+to
#TO DO separate preparation of R to avoid having it twice when from+to
