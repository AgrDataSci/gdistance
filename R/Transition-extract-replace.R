#' Extracting and replacing: class Transition
#'
#' Methods for functions \code{[} and \code{[<-} for object of
#' the class TransitionLayer.
#'
#' Methods for functions \code{[[} and \code{[[<-} for object of
#' the class TransitionStack.
#'
#' Also see \code{\link[gdistance]{adjacencyFromTransition}}.
#'
#' @name Transition-extract-replace
#' @rdname Transition-extract-replace
#' @docType methods
#' @keywords methods
#' @aliases [,TransitionLayer,numeric,numeric,missing-method
#' @aliases [,TransitionLayer,matrix,missing,missing-method
#' @aliases [<-,TransitionLayer,matrix,missing,ANY-method
#' @aliases [<-,TransitionLayer,numeric,numeric,ANY-method
#' @aliases [[,TransitionStack,numeric,missing-method
#' @aliases [[<-,TransitionStack,numeric,missing,TransitionData-method
#' @aliases TransitionStack<-
#' @aliases transitionMatrix<-
#' @aliases transitionMatrix<-,TransitionLayer,sparseMatrix-method
#'
#' @param x an object of class \code{Transition*}
#' @param value the value to assign
#' @examples
#' #Create a new raster and set all its values to unity.
#' r <- raster(nrows=18, ncols=36)
#' r <- setValues(r,rep(1,ncell(r)))
#'
#' #Create TransitionLayer objects
#' tr1 <- transition(r,mean,4)
#' tr2 <- tr1
#'
#' #Extracting and replacing
#' tr1[cbind(1:9,1:9)] <- tr2[cbind(1:9,1:9)]
#' tr1[1:9,1:9] <- tr2[1:9,1:9]
#' tr1[1:5,1:5]
#' @exportMethod transitionMatrix<-
#' @exportMethod [<-
#' @exportMethod [
#' @exportMethod [[<-
#' @exportMethod [[
setGeneric("transitionMatrix<-",
           function(x, value){
             standardGeneric("transitionMatrix<-")
             })

setReplaceMethod("transitionMatrix",
                 signature(x = "TransitionLayer",
                           value = "sparseMatrix"),
	function(x, value)
	{
		if(dim(value)[1] != dim(value)[2]){stop("sparse matrix has to be square")}
		if(dim(value)[1] != ncell(x)[2]){
		  stop("sparse matrix has to have ncell(x) rows and columns")
		}
		x@transitionMatrix <- value
		x@transitionCells <- 1:ncell(x)
		return(x)
	}
)

setGeneric("transitionMatrix<-", function(x, value) {
  standardGeneric("transitionMatrix<-")
})

setReplaceMethod ("transitionMatrix",
                  signature(x = "TransitionLayer",
                            value = "sparseMatrix"),
                  function(x, value){
                    if(dim(value)[1] != dim(value)[2]){
                      stop("sparse matrix has to be square")
                    }

                    if(dim(value)[1] == ncell(x)){x@transitionMatrix <- value}
                    else
                    {
                      if(dim(value)[1] == length(transitionCells(x)))
                      {
                        trC <- transitionCells(x)
                        tr <- Matrix(0,ncell(x),ncell(x))
                        tr[trC,trC] <- value
                        x@transitionMatrix <- tr
                      }
                      else{stop("value is of wrong dimensions; either ncell(transition)",
                                " or length(transitionCells(transition))")}
                    }
                    return(x)
                  }
)

setMethod("[", signature(x = "TransitionLayer", i="numeric",
                         j="numeric", drop="missing"), function(x,i,j)
                         {
                           tm <- transitionMatrix(x)
                           tm <- tm[i,j]
                           return(tm)
                         }
)

setMethod("[", signature(x = "TransitionLayer", i="matrix",
                         j="missing", drop="missing"), function(x,i)
                         {
                           tm <- transitionMatrix(x)
                           tm <- tm[i]
                           return(tm)
                         }
)

setMethod("[<-", signature(x = "TransitionLayer",
                           i="matrix", j="missing", value="ANY"),
          function(x, i, value){
            tm <- transitionMatrix(x)
            tm[i] <- value
            x@transitionMatrix <- tm
            return(x)
          }
)

setMethod("[<-", signature(x = "TransitionLayer",
                           i="numeric", j="numeric", value="ANY"),
          function(x, i, j, value)
          {
            tm <- transitionMatrix(x)
            tm[i,j] <- value
            transitionMatrix(x) <- tm
            return(x)
          }
)


setMethod("[[", signature(x = "TransitionStack",
                          i="numeric", j="missing"),
          function(x,i) {
  if (!(all(i %in% 1:nlayers(x)))){stop("indices should correspond to layers")}
  else
  {
    if(length(i)==1)
    {
      result <- new("TransitionLayer", nrows=as.integer(nrow(x)),
                    ncols = as.integer(ncol(x)),
                    extent = extent(c(xmin(x), xmax(x),
                                      ymin(x), ymax(x))),
                    crs=projection(x, asText=FALSE))
      result@transitionMatrix <- x@transition[[i]]@transitionMatrix
      result@transitionCells <- x@transition[[i]]@transitionCells
      result@matrixValues <- x@transition[[i]]@matrixValues
    }
    if(length(i)>1)
    {
      result <- x
      result@transition <- x@transition[i]
    }
  }
  return(result)
}
)

setMethod("[[<-", signature(x = "TransitionStack", i="numeric",
                            j="missing", value="TransitionData"),
          function(x,i, value) {
            x@transition[[i]] <- value
            return(x)
          }
)

