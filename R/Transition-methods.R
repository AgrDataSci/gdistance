#' Arithmetic and mathematical operations with objects of Transition* classes
#' 
#' Standard arithmetic operators for computations with \code{Transition*} 
#' objects and numeric values. Transition objects must have the same extent
#' and resolution. All arithmetic and mathematical operations that work on 
#' the sparse matrices are available for \code{Transition*} objects. 
#' 
#' @name ArithMath-methods
#' @docType methods
#' @keywords methods
#' @keywords math
#' @keywords spatial
#' 
#' @aliases Arith-methods
#' @aliases Arith,ANY,TransitionLayer-method
#' @aliases Arith,TransitionLayer,ANY-method
#' @aliases Arith,TransitionLayer,TransitionLayer-method
#' @aliases Arith,TransitionLayer,TransitionStack-method
#' @aliases Arith,TransitionStack,TransitionLayer-method
#' @aliases Arith,TransitionStack,ANY-method
#' @aliases Arith,ANY,TransitionStack-method
#' @aliases Arith,TransitionLayer,TransitionStack-method
#' @aliases Math-methods
#' @aliases Math,TransitionLayer-method
#' @aliases Math,TransitionStack-method
#' 
#' @param e1 objects
#' @param e2 objects
#' @return  \code{Transition*} object or numeric.
#' @examples 
#' #create a new raster and set all its values to unity.
#' raster <- raster(nrows=18, ncols=36) 
#' raster <- setValues(raster,rep(1,ncell(raster)))
#' 
#' #create TransitionLayer objects
#' tr1 <- transition(raster,mean,4)
#' tr2 <- tr1
#' 
#' #arithmetic operations
#' tr3 <- tr1 * tr2
#' tr4 <- tr3 * 4
#' 
#' #mathematical operations
#' tr5 <- sqrt(tr4)
#' @author Jacob van Etten
#' @exportMethod [
#' @exportMethod [<-
#' @exportMethod ==
#' @exportMethod coerce
#' @exportMethod Ops
#' @exportMethod Compare
#' @exportMethod Logic
#' @exportMethod Arith
#' @exportPattern "^[^\\.]"
setMethod("Arith", signature(e1 = "TransitionLayer", e2 = "TransitionLayer"),
		function(e1, e2)
		{
			if(as(e1, "BasicRaster") == as(e2, "BasicRaster"))
				{
					matrix.dsC <- callGeneric(transitionMatrix(e1), 
					                          transitionMatrix(e2))
					transitionMatrix(e1) <- matrix.dsC
					e1@transitionCells <- 1:ncell(e1)
					return(e1)
				}
			else {stop("objects do not coincide in resolution, extent and/or projection")}
		}
)

setMethod("Arith", signature(e1 = "TransitionLayer", e2 = "ANY"),
		function(e1, e2)
		{
			matrix.dsC <- callGeneric(transitionMatrix(e1,inflate=FALSE),e2)
			transitionMatrix(e1) <- matrix.dsC
			return(e1)
		}
)

setMethod("Arith", signature(e1 = "ANY", e2 = "TransitionLayer"),
		function(e1, e2)
		{
			matrix.dsC <- callGeneric(e1,transitionMatrix(e2,inflate=FALSE))
			transitionMatrix(e1) <- matrix.dsC
			return(e1)
		}
)

setMethod("Math", signature(x = "TransitionLayer"),
		function(x)
		{
			transitionMatrix(x) <- callGeneric(transitionMatrix(x, inflate=FALSE))
			return(x)
		}
)

setMethod("==", signature(e1 = "TransitionLayer", e2 = "TransitionLayer"),
		function(e1, e2)
		{
			c1 <- transitionMatrix(e1) == transitionMatrix(e1)
			c2 <- as(e1, "BasicRaster") == as(e2, "BasicRaster")
			cond <- c1 & c2
			return(cond)
		}
)

#TransitionStack
setMethod("Arith", signature(e1 = "TransitionLayer", e2 = "TransitionStack"),
		function(e1, e2)
		{
			if(as(e1, "BasicRaster") == as(e2, "BasicRaster"))
				{
					for(i in 1:nlayers(e2))
					{
						matrix.dsC <- callGeneric(as(e1,"sparseMatrix"),as(e2@transition[[i]],"sparseMatrix"))
						e2@transition[[i]]@transitionMatrix <- matrix.dsC
						e2@transition[[i]]@transitionCells <- 1:ncell(e2)
					}
					return(e2)
				}
			else {stop("objects do not coincide in resolution, extent and/or projection")}
		}
)

setMethod("Arith", signature(e1 = "TransitionStack", e2 = "TransitionLayer"),
		function(e1, e2)
		{
			result <- callGeneric(e2,e1)
			return(result)
		}
)

setMethod("Arith", signature(e1 = "TransitionStack", e2 = "ANY"),
		function(e1, e2)
		{
			if(length(e2) == 1)
			{
				for(i in 1:nlayers(e1))
				{
					matrix.dsC <- callGeneric(transitionMatrix(e1@transition[[i]]),e2)
					e1@transition[[i]]@transitionMatrix <- matrix.dsC
					e1@transition[[i]]@transitionCells <- 1:ncell(e1)
				}
				
			}
			else if(length(e2) == nlayers(e1))
			{
				for(i in 1:nlayers(e1))
				{
					matrix.dsC <- callGeneric(transitionMatrix(e1@transition[[i]]),e2[i])
					e1@transition[[i]]@transitionMatrix <- matrix.dsC
					e1@transition[[i]]@transitionCells <- 1:ncell(e1)
				}
			}
			else
			{
				stop("length of argument should be either 1 or nlayers(TransitionStack)")
			}
			return(e1)
		}
)

setMethod("Arith", signature(e1 = "ANY", e2 = "TransitionStack"),
		function(e1, e2)
		{
			result <- callGeneric(e2,e1)
			return(result)
		}
)

setMethod("Math", signature(x = "TransitionStack"),
		function(x)
		{
			for(i in 1:nlayers(e1))
			{
				matrix.dsC <- callGeneric(transitionMatrix(x@transition[[i]]))
				e1@transition[[i]]@transitionMatrix <- matrix.dsC
				e1@transition[[i]]@transitionCells <- 1:ncell(e1)
			}
			return(x)
		}
)

setMethod("==", signature(e1 = "TransitionStack", e2 = "TransitionStack"),
		function(e1, e2)
		{
			c1 <- e1@transition == e2@transition
			c2 <- as(e1, "BasicRaster") == as(e2, "BasicRaster")
			cond <- c1 & c2
			return(cond)
		}
)

#TO DO matrixValues == "resistance"
