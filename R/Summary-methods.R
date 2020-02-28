#' Summary methods
#' 
#' The following summary methods are available: mean, 
#'  Median, max, min, range, prod, sum, any, all
#' 
#' @name Summary-methods
#' @aliases Summary-methods
#' @aliases Summary,TransitionLayer-method
#' @aliases Summary,TransitionStack-method
#' @aliases sum,TransitionStack-method
#' @aliases mean,TransitionStack-method
#' @aliases summary,TransitionLayer-method
#' @aliases summary,TransitionStack-method
#' 
#' @param x objects
#' @param ... further arguments passes to or from methods
#' @param na.rm logical, should missing values be removed?
#' 
#' @return a TransitionLayer
#' 
#' @note
#'  These methods compute a summary statistic based on cell values of 
#'  layers in a TransitionStack. The result of these methods is always a 
#'  single TransitionLayer. 
#' @author Jacob van Etten
#' @examples 
#' #Create a new raster and set all its values to unity.
#' raster <- raster(nrows=18, ncols=36) 
#' raster <- setValues(raster,rep(1,ncell(raster)))
#' 
#' #Create a Transition object from the raster
#' tr <- transition(raster,mean,4)
#' 
#' trS <- stack(tr, tr*2)
#' 
#' #Apply a Summary method
#' trSum <- sum(trS)
#' 
#' #plot(raster(trMean))
#' @exportMethod Summary
setMethod("Summary", signature(x='TransitionStack'),
	function(x, ..., na.rm=FALSE){
		objectList <- list(x, ...)
		if(!is.null(objectList))
		{
			if(length(objectList)>1) {
			  warning("operations with more than two", 
			          " Transition* objects not implemented; use stack()")
			  }
		}
		call <- sys.call()
        fun <- as.character(call[[1L]])
		result <- x[[1]] * 0
		trResult <- .rowWiseFun(x,result,fun)
		transitionMatrix(result) <- trResult
		return(result)
	}
)

.rowWiseFun <- function(x, result, fun)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nrow(trResult))
		{
			rows <- x[[1]]@transitionMatrix[i,]
			for(j in 2:nlayers(x)) {rows <- rbind(rows, x[[j]]@transitionMatrix[i,])}
			trResult[i,] <- apply(rows, 2, fun)
		}
		return(trResult)
}

setMethod("Summary", signature(x = "TransitionLayer"),
		function(x, ..., na.rm)
		{
			objectList <- list(x, ...)
			if(!is.null(objectList))
			{
				if(length(objectList)>1)
				{
					stop("operations with more than two Transition*", 
					     " objects not implemented; use stack()")
				}
				else{result <- callGeneric(x@transitionMatrix,
				                           objectList[[1]]@transitionMatrix,
				                           na.rm=na.rm)}
			}
			else{result <- callGeneric(x@transitionMatrix, na.rm=na.rm)}

			return(result)
		}
)

setMethod("sum", signature(x='TransitionStack'),
	function(x, ..., na.rm=FALSE){
		objectList <- list(x, ...)
		if(!is.null(objectList))
		{
			if(length(objectList)>1){
			  warning("operations with more than two Transition*",
			          " objects not implemented; use stack()")
			  }
		}
		result <- x[[1]] * 0
		trResult <- .MatrixSum(x,result) 
		transitionMatrix(result) <- trResult
		return(result)
	}
)

.rowWiseSum <- function(x, result)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nrow(trResult))
		{
			rows <- x[[1]]@transitionMatrix[i,]
			for(j in 2:nlayers(x)) {
			  rows <- rbind(rows, x[[j]]@transitionMatrix[i,])
			}
			
			trResult[i,] <- colSums(rows)
		}
		return(trResult)
}

.MatrixSum <- function(x, result)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nlayers(x))
		{
			m <- x[[i]]@transitionMatrix
			trResult <- trResult + m
		}
		return(trResult)
}

setMethod("mean", signature(x='TransitionStack'),
	function(x, ..., na.rm=FALSE){
		objectList <- list(x, ...)
		if(!is.null(objectList))
		{
			if(length(objectList)>1){ 
			  warning("operations with more than two Transition*",
			          " objects not implemented; use stack()")
			  }
		}
		result <- x[[1]] * 0
		trResult <- .rowWiseMean(x,result)
		transitionMatrix(result) <- trResult
		return(result)
	}
)

.rowWiseMean <- function(x, result)
{
		trResult <- transitionMatrix(result)
		for(i in 1: nrow(trResult))
		{
			rows <- x[[1]]@transitionMatrix[i,]
			for(j in 2:nlayers(x)) {rows <- rbind(rows, 
			                                      x[[j]]@transitionMatrix[i,])}
			trResult[i,] <- colMeans(rows)
		}
		return(trResult)
}

