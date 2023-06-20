#' Coords class
#'
#' This is a class union, providing a convenient class for
#'  coordinates in several formats.
#'
#' The class accepts coordinates in any of the following formats:
#'  1. SpatialPoints
#'  2. two-columned matrix
#'  3. vector of length 2
#'
#' @name Coords class
#' @aliases Coords-class
#' @examples
#' showClass("Coords")
#' @keywords classes
#' @exportClass Coords
setClassUnion("Coords",
              c("numeric", "matrix", "SpatialPoints"))

