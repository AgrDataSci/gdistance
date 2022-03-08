#' RasterLayer from TransitionLayer object
#' 
#' Create a RasterLayer from a TransitionLayer with a call to 
#' the generic function \code{raster}. The n x n transition matrix of 
#' the TransitionLayer is transformed to form the values n cells of a raster.
#' @rdname rasterFromTransition
#' @docType methods
#' @aliases raster,TransitionLayer-method
#' @param x an object of class \code{Transition*}
#' @param reduceMethod character for the method to reduce the transition matrix.
#'  See details
#' @details 
#' The following methods to \sQuote{reduce} the transition matrix are
#' available with the optional argument \code{reduceMethod}):
#' \itemize{
#'     \item{colSums}
#'     \item{rowSums}
#'     \item{colMeans}
#'     \item{rowMeans}
#'     \item{NZcolMeans}
#'     \item{NZrowMeans}
#' }
#' The latter two methods only take into account the non-zero entries 
#' in the transition matrix.
#' The default is NZcolMeans.
#' @return a RasterLayer
#' @keywords spatial
#' @author Jacob van Etten
#' @examples 
#' #create a new raster and set all its values to unity.
#' r <- raster(nrows=18, ncols=36) 
#' r <- setValues(r,runif(ncell(r),0,1))
#' 
#' #create a Transition object from the raster
#' tr1 <- transition(r,mean,8)
#' 
#' #asymmetric
#' asf <- function(x) max(x) - x[1] + x[2]
#' tr2 <- transition(r,asf,8, symm=FALSE)
#' 
#' #create RasterLayer objects
#' r1 <- raster(tr1)
#' r2 <- raster(tr2)
#' r3 <- raster(tr1, "colMeans")
#' @exportMethod raster
setMethod('raster', signature(x='TransitionLayer'), 
		function(x, reduceMethod="NZcolMeans") {
		rs <- as(x,"RasterLayer")
		dataVector <- vector(length=ncell(x))
		m <- as(x,"sparseMatrix")
		if(reduceMethod == "colSums") dataVector <- colSums(m)
		if(reduceMethod == "rowSums") dataVector <- rowSums(m)
		if(reduceMethod == "colMeans") dataVector <- colMeans(m)
		if(reduceMethod == "rowMeans") dataVector <- rowMeans(m)
		if(reduceMethod == "NZcolMeans" | reduceMethod == "NZrowMeans"){
			mL <- as(m,"lMatrix")
			if(reduceMethod == "NZrowMeans") dataVector <- rowSums(m)/rowSums(mL)
			if(reduceMethod == "NZcolMeans") dataVector <- colSums(m)/colSums(mL)
		}
		rs <- setValues(rs, dataVector) 
		return(rs)
	}
)

