#' Coordinates to Matrix
#'
#' @param Coords an object with coordinates
#' @return an object of class matrix
#' @examples
#' # Write example
#' @noRd
.coordsToMatrix <- function(Coords)
{
  #  if (class(Coords) == "numeric") {
	if (!is.matrix(Coords) & is.numeric(Coords)) {
		if (length(Coords) == 2) {
		  Coords <- t(as.matrix(Coords))
		} 
		else{stop("coordinates given as a vector, 
		          but the vector does not have a length of two")}
	}
	#if (class(Coords) == "matrix") {
	if(is.matrix(Coords)) {
	  #if (!(ncol(Coords) == 2)) {
		if(dim(Coords)[[2]] != 2) {
		stop("coordinates given as a matrix,", 
		     " but the matrix does not have two columns")
		}
	}	

	if(inherits(Coords, "SpatialPoints")) {
		Coords <- sp::coordinates(Coords)
	}
	return(Coords)
}

#' Connected components
#'
#' @param x an object with coordinates
#' @return the connected components
#' @examples
#' # Write example
#' @noRd
.connected.components <- function(x)
{
	adj.graph <- igraph::graph.adjacency(transitionMatrix(x))
	
	clustermembership <- cbind(1:ncell(x), 
	                           as.integer(igraph::clusters(adj.graph)$membership) + 1)
	
	return(clustermembership)
}

#' Current components
#'
#' @param L a L object
#' @param Lr a Lr object
#' @param A an A object
#' @param n a n object
#' @param indexFrom the class 
#' @param indexTo the class to coerce
#' @return current components
#' @examples
#' # Write example
#' @noRd
.current <- function(L, Lr, A, n, indexFrom, indexTo) 
{
  # This should avoid too big floating points as "Voltage differences"
  C <- 1e-300 * n 
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- C
 	e[indexTo,] <- -C
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	#I = V * Conductance
	Current <- colSums(V * -L)/2 
	Current[indexFrom] <- 1
	Current[indexTo] <- 1
	return(Current)
}

#' Current R components
#'
#' @param L a L object
#' @param Lr a Lr object
#' @param A an A object
#' @param n a n object
#' @param indexFrom the class 
#' @param indexTo the class to coerce
#' @return current components
#' @examples
#' # Write example
#' @noRd
.currentR <- function(L, Lr, A, n, indexFrom, indexTo)
{
	lf <- length(indexFrom)
	lt <- length(indexTo)
	C <- 1e-300 * n
	# This should avoid too big floating points as "Voltage differences"
	Cf <- C / lf 
	Ct <- C / lt
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- Cf
 	e[indexTo,] <- -Ct
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	#I = V * Conductance
	Current <- colSums(V * -L) / 2 
	Current[indexFrom] <- 1
	Current[indexTo] <- 1
	return(Current)
}

#' Potential
#'
#' @param L a L object
#' @param Lr a Lr object
#' @param A an A object
#' @param n a n object
#' @param indexFrom the class 
#' @param indexTo the class to coerce
#' @return potential
#' @examples
#' # Write example
#' @noRd
.potential <- function(L, Lr, A, n, indexFrom, indexTo) 
{
  # This should avoid too big floating points as "Voltage differences"
  C <- 1e-300 * n 
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- C
 	e[indexTo,] <- -C
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	return(V)
}

#' Current M
#'
#' @param L a L object
#' @param Lr a Lr object
#' @param A an A object
#' @param n a n object
#' @param indexFrom the class 
#' @param indexTo the class to coerce
#' @param index the index
#' @return current M
#' @examples
#' # Write example
#' @noRd
.currentM <- function(L, Lr, A, n, indexFrom, indexTo, index) 
{
  # This should avoid too big floating points as "Voltage differences"
  C <- 1e-300 * n 
	e <- matrix(0, ncol=1, nrow=n)
	e[indexFrom,] <- C
 	e[indexTo,] <- -C
	x <- solve(Lr,e)
	x <- as.vector(x)
	Lplusallrows <- c(x,x[length(x)]) / C
	V1 <- A * Lplusallrows
	V2 <- t(t(A) * Lplusallrows)
	V <- abs(V1 - V2)
	# I = V * Conductance
	Current <- V[index] * -L[index] 
	return(Current)
}

#' Laplacian differential operator
#'
#' @param x a x object
#' @return Laplacian differential operator
#' @examples
#' # Write example
#' @noRd
.Laplacian <- function(x) 
{
	Laplacian <- Matrix::Diagonal(x = colSums(
	  transitionMatrix(x, inflate=FALSE))
	  ) - transitionMatrix(x, inflate=FALSE)
	
	Laplacian <- methods::as(Laplacian, "symmetricMatrix")
	
	return(Laplacian)
}

#' Transition Solidify
#'
#' @param x a x object
#' @return transition
#' @examples
#' # Write example
#' @noRd
.transitionSolidify <- function(x) {
	selection <- which(rowMeans(transitionMatrix(x, inflate=FALSE)) > 1e-300)
	x@transitionCells <- x@transitionCells[selection]
	x@transitionMatrix <- transitionMatrix(x,inflate=FALSE)[selection,selection]
	return(x)
}

#' Determine place in dist vector given place in dist matrix 
#' -- from gdistanalyst
#'
#' @param i a i object
#' @param j a j object
#' @param n a n object
#' @return distance index
#' @examples
#' # Write example
#' @noRd
.distIndex <- function(i,j,n){
  
  n*(j-1) - j*(j-1)/2 + i-j

}

#' Determine place in dist matrix given place in dist vector
#' 
#' from gdistanalyst -- should be possible speed up!
#'
#' @param i a i object
#' @param n a n object
#' @return distance matrix
#' @examples
#' # Write example
#' @noRd
.matrIndex <- function(i, n) {
	cc <- cumsum(seq((n-1),1))
	out <- matrix(nrow=length(i),ncol=2)
	for(index in 1:length(i))
	{
		out[index,2] <- min(which((cc-i[index])>=0))
		out[index,1] <- -c(0,cc)[out[index,2]] + i[index] + out[index,2]
	}
	return(out)
}


