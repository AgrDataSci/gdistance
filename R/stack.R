#' @aliases stack,TransitionLayer-method
#' @aliases stack,TransitionStack-method
#' @author Jacob van Etten
setMethod("stack", signature(x="TransitionLayer"), function(x, ...) 
{
  newStack <- as(x, "TransitionStack")
  objectList <- list(...)
  TData <- .createTData(x, objectList)
  newStack@transition <- TData
  return(newStack)	
} 
)

setMethod("stack", signature(x='TransitionStack'), function(x, ...) 
{
  newStack <- as(x, "TransitionStack")
  objectList <- list(...)
  TData <- .createTData(x, objectList)
  newStack@transition <- TData
  return(newStack)	
} 
)

.createTData <- function(x, objectList)
{
  nobj <- length(objectList)
  if(nobj<1) {stop("more than one object is needed to stack")}
  TD <- transitionData(x)
  for(i in 1:nobj) {TD <- c(TD,transitionData(objectList[[i]]))}
  return(TD)
}
