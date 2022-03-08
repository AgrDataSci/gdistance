#' Reciprocal of the sum of the reciprocals of conductance values in Transition* objects
#' 
#' Reciprocal of the sum of the reciprocals of conductance Transition* objects
#' 
#' @name sumReciprocal
#' @aliases sumReciprocal
#' @keywords spatial
#' 
#' @param x1 \code{TransitionLayer} object
#' @param x2 \code{TransitionLayer} object
#' @return \code{TransitionLayer} object containing conductance values.
#' @details 
#' To calculate the total resistance of two resistors that are serially 
#'  connected, we should add their resistance values. However, if we work
#'  with conductance values, we need to take the reciprocal of the summed 
#'  reciprocals of the conductance values. This function does that when 
#'  adding two TransitionLayers with conductance values 
#'  (\code{matrixValues(tr) == "conductance"}). 
#' 
#' For a TransitionLayer with resistance values 
#'  (\code{matrixValues(tr) == "resistance"}), the function will not take
#'  reciprocals for that object, but will still take a reciprocal for the 
#'  final product (which will consequently have conductance values).
#' @examples 
#' #Create a new raster and set all its values to unity.
#' raster <- raster(nrows=18, ncols=36) 
#' raster <- setValues(raster,rep(1,ncell(raster)))
#' 
#' #Create TransitionLayer objects
#' tr1 <- transition(raster,mean,4)
#' tr2 <- tr1
#' matrixValues(tr1)
#' 
#' #Set one to resistance
#' matrixValues(tr2) <- "resistance"
#' 
#' #Sum the two objects
#' sumReciprocal(tr1,tr2)
#' @export
sumReciprocal <- function(x1, x2)
{
	if(matrixValues(x1) == "conductance") 
	{
		x1@transitionMatrix@x <- 1 / x1@transitionMatrix@x
		x1@transitionMatrix@x[x1@transitionMatrix@x == Inf] <- 0
	}
	if(matrixValues(x2) == "conductance") 
	{
		x2@transitionMatrix@x <- 1 / x2@transitionMatrix@x
		x2@transitionMatrix@x[x2@transitionMatrix@x == Inf] <- 0
	}
	newTransition <- x1 + x2
	newTransition@transitionMatrix@x <- 1 / newTransition@transitionMatrix@x
	if(sum(newTransition@transitionMatrix@x == Inf) > 0) 
	{
		warning("Inf values introduced. Set to 0.")
		newTransition@transitionMatrix@x[newTransition@transitionMatrix@x == Inf] <- 0
	}
	matrixValues(newTransition) <- "conductance"
	return(newTransition)
}

