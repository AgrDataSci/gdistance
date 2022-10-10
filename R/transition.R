#' Create an object of the class Transition
#' 
#' Create a Transition object from a RasterLayer or RasterBrick object.
#'  Transition values are calculated with a user-defined function from 
#'  the grid values.
#'
#' @name transition
#' @aliases transition
#' @aliases transition,RasterLayer-method
#' @aliases transition,RasterBrick-method
#' @keywords spatial
#' @family Transition
#' 
#' @param x \code{RasterLayer} or \code{RasterBrick} (raster package)
#' @param transitionFunction Function to calculate transition values 
#'  from grid values
#' @param directions Directions in which cells are connected 
#'  (4, 8, 16, or other), see \code{\link[raster:adjacent]{adjacent}}
#' @param ... additional arguments, passed to methods
#' 
#' @details 
#' 
#' Users may use one of three methods to construct a Transition* 
#' object with this function.
#' 
#' 1) \code{TransitionLayer} from \code{RasterLayer}
#'
#' \code{transition(x, transisitonFunction, directions, symm)}
#' 
#' When a symmetric transition matrix is required, the user should
#' supply a transitionFunction f that obeys f(i,j) = f(j,i) 
#' (a commutative function). 
#' 
#' The function \code{transition} does no commutativity check.
#' 
#' To obtain an asymmetric transition matrix, a non-commutative 
#' function should be supplied and an additional argument `symm' 
#' should be set to FALSE.
#' 
#' 2) \code{TransitionLayer} from \code{RasterBrick}
#' 
#' \code{transition(x, transitionFunction="mahal", directions)}
#' 
#' This method serves to summarize several layers of data in a single 
#' distance measure. The distance between adjacent cells is the normalized 
#' reciprocal of the Mahalanobis distance 
#' (mean distance / (mean distance + distance ij).
#' 
#' 3) \code{TransitionStack} from \code{RasterLayer}
#' 
#' In contrast with the above methods, this method produces resistance 
#' matrices by default.
#' 
#' a) Continuous variables - barriers
#' 
#' \code{transition(x, transitionFunction="barriers", 
#' directions, symm, intervalBreaks)}
#' 
#' This method creates a \code{TransitionStack} with each layer 
#' containing a discrete boundary between areas in \code{x}. 
#' Areas are defined by intervals in \code{x}.
#' The argument \code{intervalBreaks} is a vector of interval 
#' breaks corresponding to the values in \code{x}.
#' If between a pair of cells i and j, min(i,j) < break AND max(i,j) > break, 
#' then the value ij in the transition matrix becomes 1.
#' 
#' All other values in the transition matrix remain 0.
#' The package classInt offers several methods to define intervals.
#' If symm is changed from the default (TRUE) to "up" or "down", 
#' it will give either only the upslope (symm="up") or 
#' downslope (symm="down") barriers.
#' 
#' b) Categorical variables - barriers 
#' 
#' \code{transition(x, transitionFunction="barriers", directions)}
#' 
#' In this case, areas are defined as categories in the input raster.
#' A raster with a categorical variable can be created with 
#' \code{asFactor()}.
#' The layers of the resulting TransitionStack contain 
#' all possible combinations of categories.
#' Which layer contains the combination of categories i and j 
#' out of n categories, can be determined with these formulae:
#' 
#' if \code{symm} is \code{TRUE}: layer(i,j) = n*(j-1) - j*(j-1)/2 + i-j.
#' if \code{symm} is \code{FALSE} and i>j: layer(i,j) = ((n*(j-1) - j*(j-1)/2 + i-j) * 2) - 1.
#' if \code{symm} is \code{FALSE} and i<j: layer(i,j) = (n*(j-1) - j*(j-1)/2 + i-j) * 2.
#' 
#' c) Categorical variables - areas
#' 
#' \code{transition(x, transitionFunction="areas", directions)}
#' 
#' Here, areas are also a categorical variable (see under 3b).
#' The layers in the resulting TransitionStack represent each one area.
#' Connections between two cells which are each inside the area are set to 1.
#' Connections between a cell inside and a cell outside the area are set to 0.5.
#' Connections between two cells outside the area are set to 0.
#' @exportMethod transition
setGeneric("transition", function(x, transitionFunction, directions, ...) standardGeneric("transition"))

#' @exportMethod transition
setMethod("transition", 
          signature(x = "RasterLayer"), 
          def = function(x, transitionFunction, directions, 
                         symm=TRUE, intervalBreaks=NULL)
		{
			if(is(transitionFunction, "character"))	{
				if(transitionFunction != "barriers" & transitionFunction != "areas") {
					stop("argument transitionFunction invalid")
				}
				if(transitionFunction=="barriers") {
					return(.barriers(x, directions, symm, intervalBreaks))
				}
				if(transitionFunction=="areas")
				{
					return(.areas(x, directions))
				}
			} else {
				return(.TfromR(x, transitionFunction, directions, symm))
			}
		}
)

.TfromR <- function(x, transitionFunction, directions, symm)
{
	tr <- new("TransitionLayer",
		nrows=as.integer(nrow(x)),
		ncols=as.integer(ncol(x)),
		extent=extent(x),
		crs=projection(x, asText=FALSE),
		transitionMatrix = Matrix(0,ncell(x),ncell(x)),
		transitionCells = 1:ncell(x))
	transitionMatr <- transitionMatrix(tr)
	Cells <- which(!is.na(getValues(x)))
	adj <- adjacent(x, cells=Cells, pairs=TRUE, 
	                target=Cells, 
	                directions=directions)
	if(symm){adj <- adj[adj[,1] < adj[,2],]}
	dataVals <- cbind(getValues(x)[adj[,1]],
	                  getValues(x)[adj[,2]])
	transition.values <- apply(dataVals,1,transitionFunction)
	
	if(!all(transition.values>=0)){
	  warning("transition function gives negative values")
	}
	
	transitionMatr[adj] <- as.vector(transition.values)
	if(symm)
	{
		transitionMatr <- forceSymmetric(transitionMatr)
	}
	transitionMatrix(tr) <- transitionMatr
	matrixValues(tr) <- "conductance"
	return(tr)
}

.barriers <- function(x, directions, symm, intervalBreaks) {
	Xlayer <- new("TransitionLayer",
		nrows=as.integer(nrow(x)),
		ncols=as.integer(ncol(x)),
		extent=extent(x),
		crs=projection(x, asText=FALSE),
		transitionMatrix = Matrix(0,ncell(x),ncell(x)),
		transitionCells = 1:ncell(x))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)
	
	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)
		
		if(symm)
		{
			maxn <- (n^2 - n)/2
			for(i in 1:maxn)
			{
				j <- .matrIndex(i,n)
				XlayerNew <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacent(x, cells=cells1, pairs=TRUE, 
				                 target=cells2, 
				                 directions=directions)
				adj2 <- adjacent(x, cells=cells2, 
				                 pairs=TRUE,
				                 target=cells1, 
				                 directions=directions)
				adj <- rbind(adj1,adj2)
				XlayerNew[adj] <- 1
				Xstack <- stack(Xstack, XlayerNew)
			}
		} else {
			maxn <- (n^2 - n)/2
			for(i in 1:maxn)
			{
				j <- .matrIndex(i,n)
				XlayerNew1 <- Xlayer
				XlayerNew2 <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacent(x, 
				                 cells=cells1,
				                 pairs=TRUE, 
				                 target=cells2, 
				                 directions=directions)
				adj2 <- adjacent(x, 
				                 cells=cells2, 
				                 pairs=TRUE, 
				                 target=cells1, 
				                 directions=directions)
				XlayerNew1[adj1] <- 1
				XlayerNew2[adj2] <- 1				
				Xstack <- stack(Xstack, XlayerNew1, XlayerNew2)
			}
		}
	
	} else {
	
		Xmin <- transition(x, min, directions)
		Xmax <- transition(x, max, directions)
		index1 <- adjacent(x, 
		                   cells=1:ncell(x), 
		                   pairs=TRUE, 
		                   target=1:ncell(x),
		                   directions=directions)
		XminVals <- Xmin[index1]
		XmaxVals <- Xmax[index1]

		if(symm == TRUE)
		{
			for(i in 1:length(intervalBreaks))
			{
				index2 <- index1[XminVals < intervalBreaks[i] & 
				                   XmaxVals > intervalBreaks[i],]
				XlayerNew <- Xlayer
				XlayerNew[index2] <- 1
				Xstack <- stack(Xstack,XlayerNew)
			}
		}
		if(symm=="up" | symm=="down"){
		  stop("not implemented yet")
		}
		
	}
	
	Xstack <- Xstack[[2:nlayers(Xstack)]]	
	return(Xstack)
}


.areas <- function(x, directions) {

	Xlayer <- new("TransitionLayer",
		nrows=as.integer(nrow(x)),
		ncols=as.integer(ncol(x)),
		extent=extent(x),
		crs=projection(x, asText=FALSE),
		transitionMatrix = Matrix(0,ncell(x),ncell(x)),
		transitionCells = 1:ncell(x))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)
	
	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)
		
			for(i in 1:n)
			{
				transitionFunction <- function(v) {return(sum(v == i) / 2)}
				XlayerNew <- .TfromR(x, transitionFunction, directions, symm=TRUE)
				Xstack <- stack(Xstack,XlayerNew)
			}
			
	} else {
		warning("not yet implemented for raster with non-factor",
		        " variables. Contact author.")
	}
	Xstack <- Xstack[[2:nlayers(Xstack)]]	
	return(Xstack)
}


#' @export
setMethod("transition", signature(x = "RasterBrick"), 
          def = function(x, transitionFunction="mahal", directions)
		{
			if (transitionFunction != "mahal") {
				stop("only Mahalanobis distance method implemented for RasterBrick \n")
			}
			xy <- cbind(1:ncell(x), getValues(x))
			
			xy <- na.omit(xy)
			
			dataCells <- xy[,1]
			
			adj <- adjacent(x, cells = dataCells, pairs = TRUE,
			                target = dataCells, directions = directions)
			
			x.minus.y <- raster::getValues(x)[adj[, 1], ] - raster::getValues(x)[adj[, 2], ]
			
			cov.inv <- solve(cov(xy[,-1]))
			
			mahaldistance <- apply(x.minus.y,1,function(x){sqrt((x%*%cov.inv)%*%x)})
			
			mahaldistance <- mean(mahaldistance)/(mahaldistance+mean(mahaldistance))
			
			transitiondsC <- new("dsCMatrix", 
					p = as.integer(rep(0,ncell(x)+1)),
					Dim = as.integer(c(ncell(x),ncell(x))),
					Dimnames = list(as.character(1:ncell(x)),as.character(1:ncell(x))))
			
			transitiondsC[adj] <- mahaldistance
			
			nr <- as.integer(nrow(x))
			nc <- as.integer(ncol(x))
			
			tr <- new("TransitionLayer",
			          transitionMatrix = transitiondsC,
			          nrows = nr,
			          ncols = nc,
			          extent = extent(x),
			          crs = projection(x, asText = FALSE), 
			          matrixValues = "conductance")
			
			return(tr)
			
		}
)




