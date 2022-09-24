#' Transition classes
#' 
#' TransitionLayer and TransitionStack (or \code{Transition*}) are the core 
#' classes of the package gdistance. They are the main input into the 
#' functions to calculate distances and routes. An object of the class
#' TransitionLayer contains two main elements:
#' 
#' a. a transition matrix with transition values between connected 
#'  cells in a raster - an object of class sparseMatrix (package Matrix); 
#' b. information on the extent, resolution and projection of the underlying 
#'  raster - an object of class Raster (package raster).
#'  
#' All slots belong to these two elements from other package, 
#'  except two additional slots:
#'  1. slot transitionCells, which is only used internally in the package;
#'  2. slot matrixValues indicates if the nonzero values of the transition 
#'  matrix contains conductance or resistance values.
#' 
#' Class TransitionStack contains one or more transition matrices.
#' 
#' Class Transition is the union of TransitionLayer and TransitionStack.
#' 
#' @name Transition-classes
#' @aliases TransitionLayer-class
#' @aliases TransitionStack-class
#' @aliases TransitionData-class
#' @aliases Transition-class
#' @aliases coerce,TransitionLayer,sparseMatrix-method
#' @aliases coerce,TransitionLayer,RasterLayer-method
#' @aliases coerce,TransitionLayer,TransitionStack-method
#' @aliases coerce,TransitionLayer,TransitionData-method
#' @aliases coerce,RasterLayer,TransitionLayer-method
#' @aliases coerce,TransitionData,sparseMatrix-method
#' @aliases show,TransitionLayer-method
#' @aliases show,TransitionStack-method
#' @aliases ==,TransitionLayer,TransitionLayer-method
#' @aliases ==,TransitionStack,TransitionStack-method
#' @keywords classes
#' @section Objects from the Class:
#' Objects can be created by calls of the form 
#'  new("Transition", nrows, ncols, xmin, xmax, ymin, ymax, projection).
#' @section  Extends: 
#' Class \code{\linkS4class{Raster}}
#' @examples 
#' showClass("TransitionLayer")
#' 
#' tr <- new("TransitionLayer", nrows=as.integer(36), ncols=as.integer(18),
#'           extent=extent(c(xmin=-180,xmax=180, ymin=-90,ymax=90)),
#'           crs=CRS("+proj=longlat +datum=WGS84"))
#' 
#' tr <- new("TransitionLayer",nrows=as.integer(36),ncols=as.integer(18),
#'           extent=extent(c(xmin=-180, xmax=180, ymin=-90,ymax=90)),
#'           crs=CRS(""))
#'           
#' @slot transitionMatrix Object of class \code{"sparseMatrix"}
#' @slot transitionCells Object of class \code{"integer"}
#' @slot matrixValues Object of class \code{"character"}
#' @slot ncols Object of class \code{"integer"}
#' @slot nrows Object of class \code{"integer"}
#' @slot crs Object of class \code{"CRS"} from \code{sp} package
#' @slot extent Object of class \code{"Extent"}
#' @slot layername Object of class \code{"vector"}
#' @export
setClass(Class="TransitionData",
         representation = representation(
           transitionMatrix = "sparseMatrix",
           transitionCells = "numeric",
           matrixValues = "character"
         ),
         prototype(
           transitionMatrix = Matrix(0,1,1),
           transitionCells = 1,
           matrixValues = "conductance"
         ),
         validity = function(object){
             cond1 <- (nrow(object@transitionMatrix) == ncol(object@transitionMatrix))
             cond2 <- (object@matrixValues == "resistance" | object@matrixValues == "conductance")
             cond3 <- length(transitionCells(object)) == object@transitionMatrix@Dim[1]
             cond <- cond1 & cond2 & cond3 
             return(cond)
         }
)

#' @export
setClass(Class="TransitionLayer",
         contains = c("BasicRaster", "TransitionData"),
         prototype = prototype(
           rotated = FALSE,
           ncols = as.integer(1),
           nrows = as.integer(1),
           layernames = c(""),
           unit=c(""),
           z = list(),
           crs = CRS(as.character(NA)),
           transitionMatrix = Matrix(0,1,1),
           transitionCells = 1,
           matrixValues = "conductance"
         ),
         validity = function(object){
           cond1 <- (nrow(object@transitionMatrix) == ncol(object@transitionMatrix)) 
           cond2 <- (object@matrixValues == "resistance" | object@matrixValues == "conductance")
           cond3 <- length(transitionCells(object)) == object@transitionMatrix@Dim[1]
           cond <- cond1 & cond2 & cond3 
           return(cond)
         }
)

#' @export
setClass("TransitionStack",
         contains = "BasicRaster",
         representation (
           transition = "list"
         ),
         prototype(
           rotated = FALSE,
           ncols= as.integer(1),
           nrows= as.integer(1),
           layernames = c(""),
           unit=c(""),
           z = list(),
           crs = CRS(as.character(NA)),
           transition = list(new("TransitionData"))
         ),
         validity = function(object) {return(TRUE)}
)

#' @export
setClassUnion("Transition", c("TransitionLayer", "TransitionStack"))
