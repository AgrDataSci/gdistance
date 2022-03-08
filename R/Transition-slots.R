#' Extract or change elements of Transition* objects
#' 
#' These functions are to be used to access slots of Transition* objects.
#' 
#' @name Transition-slots
#' @rdname Transition-slots
#' @docType methods
#' @keywords methods
#' @aliases transitionMatrix
#' @aliases nlayers,TransitionLayer-method
#' @aliases nlayers,TransitionStack-method
#' @aliases transitionMatrix,TransitionLayer,logical-method
#' @aliases transitionMatrix,TransitionLayer,missing-method
#' @aliases transitionMatrix,TransitionData,missing-method
#' @aliases matrixValues
#' @aliases matrixValues<-
#' @aliases transitionCells
#' @aliases transitionData
#' @aliases transitionCells,TransitionLayer-method
#' @aliases transitionCells,TransitionData-method
#' @aliases matrixValues,TransitionLayer-method
#' @aliases matrixValues<-,TransitionLayer,character-method
#' @aliases matrixValues,TransitionStack-method
#' @aliases matrixValues<-,TransitionStack,character-method
#' @aliases transitionData,TransitionLayer-method
#' @aliases transitionData,TransitionStack-method
#' 
#' @param x object of class \code{Transition*}
#' @param inflate logical (default is \code{TRUE})
#' 
#' @exportMethod transitionMatrix
#' @exportMethod transitionCells
#' @exportMethod matrixValues
#' @exportMethod matrixValues<-
#' @exportMethod nlayers
#' @exportMethod transitionData
setGeneric("transitionMatrix", function(x, inflate) {
  standardGeneric("transitionMatrix")
})

setMethod("transitionMatrix",
          signature(x = "TransitionLayer", 
                    inflate="missing"),
          function(x)
          {
            transitionMatrix(x=x, inflate=TRUE)
          }
)

setMethod("transitionMatrix", 
          signature(x = "TransitionLayer",
                    inflate="logical"),
          function(x, inflate)
          {
            .tr(x, inflate, ncell(x))
          }
)

setMethod("transitionMatrix", 
          signature(x = "TransitionData", 
                    inflate="missing"),
          function(x)
          {
            .tr(x=x, inflate=FALSE, nc=0)
          }
)

.tr <- function(x, inflate, nc)
{
  if(inflate & length(transitionCells(x)) != ncell(x))
  {
    tr <- Matrix(0, nc,nc)
    cells <- transitionCells(x)
    tr[cells,cells] <- x@transitionMatrix
  }
  if(!inflate | length(transitionCells(x)) == nc)
  {
    tr <- x@transitionMatrix
  }
  return(tr)
}



setGeneric("transitionCells", function(x){
  standardGeneric("transitionCells")
}
)

setMethod ("transitionCells", signature(x = "TransitionLayer"),
           function(x)
           {
             return(x@transitionCells)
           }
)

setMethod ("transitionCells", signature(x = "TransitionData"),
           function(x)
           {
             return(x@transitionCells)
           }
)

setGeneric("matrixValues", function(x) {
  standardGeneric("matrixValues")
  })

setMethod("matrixValues", 
           signature(x = "TransitionLayer"),
           function(x){x@matrixValues}
)

setMethod("matrixValues", 
           signature(x = "TransitionStack"),
           function(x){stop("not implemented yet")}
)

setGeneric("matrixValues<-", function(x, value) {
  standardGeneric("matrixValues<-")
})

setReplaceMethod ("matrixValues", 
                  signature(x = "TransitionLayer", value = "character"),
                  function(x, value){
                    if (value == "resistance" | value == "conductance") 
                    {
                      x@matrixValues <- value
                      return(x)
                    }
                    else {stop("matrixValues can only be set to resistance or conductance")}
                  }
)

setMethod("nlayers", signature(x="TransitionStack"), 
          function(x)
          {
            return(length(x@transition)) 
          }
)

setGeneric("transitionData", function(x) {
  standardGeneric("transitionData")
})

setMethod("transitionData", signature(x = "TransitionLayer"),
           function(x){
             as(x, "TransitionData")
})

setMethod("transitionData", signature(x = "TransitionStack"),
           function(x){
             x@transition
})

