#' Randomized shortest path distance
#' 
#' Calculates the randomized shortest path distance between points.
#' 
#' @name rSPDistance
#' @aliases rSPDistance
#' @aliases rSPDistance,TransitionLayer,Coords-method
#' @keywords spatial
#' 
#' @param x \code{TransitionLayer} object
#' @param from point locations coordinates (of SpatialPoints, 
#' matrix or numeric class)
#' @param to point locations coordinates (of SpatialPoints, 
#' matrix or numeric class)
#' @param theta theta is the degree from which the path randomly
#' deviates from the shortest path, 0 < theta < 20
#' @param totalNet total or net movements between cells
#' @param method method 1 (as defined in Saerens et al.) or 
#' method 2 (a modified version, see below in Details)
#' 
#' @return distance matrix (S3 class dist or matrix)
#' @details 
#' The function implements the algorithm given by Saerens et al. (2009). 
#' 
#' Method 1 implements the method as it is. 
#' Method 2 uses W = exp(-theta * ln(P)).
#' 
#' @author Jacob van Etten
#' 
#' @references 
#' Saerens M., L. Yen, F. Fouss, and Y. Achbany. 2009. 
#' Randomized shortest-path problems: two related models.
#' \emph{Neural Computation}, 21(8):2363-2404.
#' @examples 
#' #Create a new raster and set all its values to unity.
#' r <- raster(nrows=18, ncols=36) 
#' r <- setValues(r,rep(1,ncell(raster)))
#' 
#' #Create a Transition object from the raster
#' tr <- transition(r,mean,4)
#' 
#' #Create two sets of coordinates
#' sP1 <- SpatialPoints(cbind(c(65,5,-65),c(55,35,-35)))
#' sP2 <- SpatialPoints(cbind(c(50,15,-40),c(80,20,-5)))
#' 
#' #Calculate the RSP distance between the points
#' rSPDistance(tr, sP1, sP2, 1)
#' 
#' @seealso \code{\link{geoCorrection}}
#' @export

rSPDistance <- function(x, from, to, theta, totalNet="net", method=1)
{
	if(theta < 0 | theta > 20 ) {stop("theta value out of range (between 0 and 20)")}
  if(method != 1 & method != 2) {stop("method should be either 1 or 2")}
	
	cellnri <- raster::cellFromXY(x, from)
	cellnrj <- raster::cellFromXY(x, to)
	transition <- .transitionSolidify(x)
	tc <- transitionCells(x)

	ci <- match(cellnri,tc)
	cj <- match(cellnrj,tc)
		
	tr <- transitionMatrix(x, inflate=FALSE)
	
	.rSPDist(tr, ci, cj, theta, totalNet, method)
}

.rSPDist <- function(tr, ci, cj, theta, totalNet, method)
{
	trR <- tr
	trR@x <- 1 / trR@x 
	nr <- dim(tr)[1] 
	Id <- Matrix::Diagonal(nr) 
	P <- .normalize(tr, "row")
  

  
  if(method == 1)
  {

    W <- trR
    #zero values are not relevant because of next step exp(-theta * trR@x)
	  W@x <- exp(-theta * trR@x)  
    W <- W * P
    
  }	else	{
    
    adj <- adjacencyFromTransition(tr)
	  W <- trR
	  #if the value is 1 then you get a natural random walk
	  W[adj] <- exp(-theta * -log(P[adj])) 

  }
   
  #if (any(rowSums(W) < 1)) warning("one or more row sums of W are < 1")
  
	D <- matrix(0, nrow=length(ci), ncol=length(cj))
	
	for(j in 1:length(cj))
	{
	  Ij <- Matrix::Diagonal(nr)
		Ij[cj[j],cj[j]] <- 0
		Wj <- Ij %*% W
		IdMinusWj <- methods::as((Id - Wj), "dgCMatrix")		
		ej <- rep(0,times=nr)
		ej[cj[j]] <- 1
		zcj <- solve(IdMinusWj, ej)
    
		for(i in 1:length(ci))
		{
			ei <- rep(0,times=nr)
			ei[ci[i]] <- 1
			zci <- solve(t(IdMinusWj),ei)
			zcij <- sum(ei*zcj)

			N <- (Matrix::Diagonal(nr, 
			                       as.vector(zci)) %*% Wj %*% Matrix::Diagonal(nr, 
			                                                       as.vector(zcj))) / zcij
      
      if(totalNet == "net")
      {
        #N is here the NET number of passages, like McRae-random walk
        N <- Matrix::skewpart(N) * 2 
        N@x[N@x<0] <- 0
      }

      # Computation of the cost dij between node i and node j
			D[i,j] <-  sum(trR * N)
			
    }
	}
	
	return(D)

}

