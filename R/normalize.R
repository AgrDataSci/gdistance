#' Normalize
#' 
#' Normalize the transition matrix.
#' 
#' @name normalize
#' @aliases normalize
#' @aliases normalize,TransitionLayer-method
#' @keywords spatial
#' @keywords methods
#' 
#' @param x object of class \code{Transition*}
#' @param ... optional argument \code{method} (see Details)
#' @return an object of class \code{TransitionLayer}
#' @details 
#' argument \code{method} passed through \code{...} 
#'  a character for the normalization method
#'  the default is 'row' users can set the optional method argument 
#'  to either "col" or "symm"
#'
#' Normalization of the weighted adjacency matrix in the Transition* object. 
#'  Matrix values are divided by their respective row-sums, column-sums, 
#'  or the product of the square-roots of both (symmetric normalization). 
#'  
#' The default  \code{method} is row-normalization. To use the other normalization 
#'  methods, users can set the optional \code{method} argument 
#'  to either "col" or "symm". 
#'  
#' For random walk calculations a symmetric matrix is needed (method = "symm").
#'  
#' @references 
#' von Luxburg, U. 2007. A tutorial on spectral clustering. 
#'  Statistics and Computing 17(4), 395-416. 
#'  \doi{https://doi.org/10.1007/s11222-007-9033-z}
#' 
#' Chung, F. 1997. Spectral Graph Theory. Conference Board of 
#'  the Mathematical Sciences, Washington.
#' @author Jacob van Etten
#' @examples
#' library("raster")
#' r <- raster(ncol=36,nrow=18)
#' r <- setValues(r,rep(1,times=ncell(r)))
#' tr <- transition(r, mean, directions=8)
#' 
#' normalize(tr)
#' 
#' normalize(tr, method="symm")
#' 
#' @exportMethod normalize
setGeneric("normalize", function(x, ...) {
  standardGeneric("normalize")
}
)

setMethod("normalize", 
          signature(x = "TransitionLayer"), 
          def = function(x, method="row")
{
  tr <- transitionMatrix(x)
  tr <- .normalize(tr, method)
  transitionMatrix(x) <- tr
  return(x)
}
)

.normalize <- function(x, method)
{
  
  if(!(method %in% c("row","col","symm"))){
    stop("invalid method argument")
  }
  
  if(method=="symm")
  {
    rs <- rowSums(x)^-.5
    cs <- colSums(x)^-.5
    
    tr <- x * rs
    tr <- t(tr)
    tr <- tr * cs
    
    tr <- t(tr)
    
    if(isSymmetric(x)) 
    {
      tr <- forceSymmetric(tr)
      tr <- as(tr, "CsparseMatrix")
    }
  }
  
  if(method=="row")
  {
    rs <- 1 / rowSums(x)
    rs[rs == Inf] <- 0
    tr <- x * rs
  }
  
  if(method=="col")
  {
    rs <- 1 / colSums(x)
    rs[rs == Inf] <- 0
    tr <- t(t(x) * rs)
  }
  
  return(tr)
}

