#' Shortest path
#' 
#' Calculates the shortest path from an origin to a goal
#' 
#' @name shortestPath
#' @aliases shortestPath
#' @aliases shortestPath,TransitionLayer,Coords,Coords-method
#' 
#' @param x \code{TransitionLayer} object
#' @param origin SpatialPoints, vector or matrix with coordinates, 
#' at the moment only the first cell is taken into account
#' @param goal SpatialPoints, vector or matrix with coordinates
#' @param ... Additional argument: output
#' @return Transition object
#' @details 
#' As an additional argument output, the desired output object can be specified: 
#'  either \dQuote{TransitionLayer} (default),  
#'  \dQuote{TransitionStack} or \dQuote{SpatialLines}.
#'  
#' If there is more than one path either (1) transition values in the 
#'  TransitionLayer get values of 1 / number of paths or (2) the SpatialLines 
#'  object will contain more than one line.
#' @author Jacob van Etten
#' @seealso \code{\link{costDistance}}, \code{\link{accCost}}
#' @examples 
#' #example equivalent to that in the documentation on r.cost/r.drain in GRASS
#' r <- raster(nrows=6, ncols=7, xmn=0, xmx=7, ymn=0, ymx=6, crs="+proj=utm +units=m")
#' 
#' r[] <- c(2, 2, 1, 1, 5, 5, 5,
#'          2, 2, 8, 8, 5, 2, 1,
#'          7, 1, 1, 8, 2, 2, 2,
#'          8, 7, 8, 8, 8, 8, 5,
#'          8, 8, 1, 1, 5, 3, 9,
#'          8, 1, 1, 2, 5, 3, 9)
#' 
#' tr <- transition(r, function(x) 1/mean(x), 8) 
#' # 1/mean: reciprocal to get permeability
#' tr <- geoCorrection(tr)
#' 
#' c1 <- c(5.5,1.5) 
#' c2 <- c(1.5,5.5)
#' 
#' #make a SpatialLines object for visualization
#' sPath1 <- shortestPath(tr, c1, c2, output="SpatialLines")
#' plot(r)
#' lines(sPath1)
#' 
#' #make a TransitionLayer for further calculations
#' sPath2 <- shortestPath(tr, c1, c2)
#' 
#' plot(raster(sPath2))
#' 
#' @exportMethod shortestPath
setGeneric("shortestPath", function(x, origin, goal, ...) {
  standardGeneric("shortestPath")
}
)

#check if Transition and RasterLayers coincide, etc.
#' @export
setMethod("shortestPath", 
          signature(x = "TransitionLayer", 
                    origin = "Coords", 
                    goal = "Coords"), 
          def = function(x, origin, goal, 
                         output="TransitionLayer")
	{
		origin <- .coordsToMatrix(origin)
		goal <- .coordsToMatrix(goal)
		return(.shortestPath(x, origin, goal, output))		
	}
)

.shortestPath <- function(x, origin, goal, output)
{
	originCells <- raster::cellFromXY(x, origin)
	goalCells <- raster::cellFromXY(x, goal)
	indexOrigin <- originCells 
	indexGoal <- goalCells 
  y <- transitionMatrix(x)
	if(isSymmetric(y)) {
	  mode <- "undirected"
	}else{
	    mode <- "directed"
	}
  
	adjacencyGraph <- igraph::graph.adjacency(y, mode=mode, weighted=TRUE)
	E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

	shortestPaths <- get.shortest.paths(adjacencyGraph,
	                                    indexOrigin, indexGoal)$vpath
	
	if(output=="TransitionLayer")
	{
		
		result <- x
		transitionMatrix(result) <- Matrix::Matrix(0, ncol=ncell(x), nrow=ncell(x))			
		for(i in 1:length(shortestPaths))
		{
			sPVector <- shortestPaths[[i]]
			adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
			adj <- rbind(adj,cbind(adj[,2], adj[,1]))
			transitionMatrix(result)[adj] <- 1/length(shortestPaths) + transitionMatrix(result)[adj]
		}
	}

	if(output=="TransitionStack")
	{
		result <- x
		transitionMatrix(result) <- Matrix::Matrix(0, ncol=ncell(x), nrow=ncell(x))			
		for(i in 1:length(shortestPaths))
		{
			resultNew <- result
			sPVector <- shortestPaths[[i]] 
			adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
			adj <- rbind(adj,cbind(adj[,2], adj[,1]))
			transitionMatrix(resultNew)[adj] <- 1/length(shortestPaths)
			result <- raster::stack(result, resultNew)
		}
		result <- result[[2:nlayers(result)]]	
	}

	if(output=="SpatialLines")
	{
		linesList <- vector(mode="list", length=length(shortestPaths))
				
		for(i in 1:length(shortestPaths))
		{
			sPVector <- shortestPaths[[i]]
			coords <- raster::xyFromCell(x, sPVector)
			linesList[[i]] <- Line(coords)
		}
		
  # Suggested by Sergei Petrov 
		LinesObject <- mapply(Lines, 
		                      slinelist = linesList,  
		                      ID = as.character(1:length(shortestPaths)),
		                      SIMPLIFY = FALSE)
		
		result <- sp::SpatialLines(LinesObject, proj4string = sp::CRS(projection(x)))
	}

	return(result)

}

