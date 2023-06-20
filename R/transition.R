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
#' \code{transition(x, transitionFunction = "mahal", directions)}
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
#' \code{transition(x, transitionFunction = "barriers",
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
#' it will give either only the upslope (symm = "up") or
#' downslope (symm = "down") barriers.
#'
#' b) Categorical variables - barriers
#'
#' \code{transition(x, transitionFunction = "barriers", directions)}
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
#' \code{transition(x, transitionFunction = "areas", directions)}
#'
#' Here, areas are also a categorical variable (see under 3b).
#' The layers in the resulting TransitionStack represent each one area.
#' Connections between two cells which are each inside the area are set to 1.
#' Connections between a cell inside and a cell outside the area are set to 0.5.
#' Connections between two cells outside the area are set to 0.
#' @exportMethod transition
setGeneric("transition", function(x, transitionFunction, directions, ...) standardGeneric("transition"))

#' @exportMethod transition
setMethod(
  "transition",
  signature(x = "RasterLayer"),
  function(x, transitionFunction, directions, symm = TRUE, intervalBreaks = NULL) {
    if (.isGlobalLonLat(x)) {
      message("The extent and CRS indicate this raster is a global lat/lon raster. ",
              "This means that transitions going off of the East or West edges ",
              "will 'wrap' to the opposite edge.")
    }

    if(is(transitionFunction, "character"))	{
      if(transitionFunction != "barriers" & transitionFunction != "areas") {
        stop("argument transitionFunction invalid")
      }
      if(transitionFunction == "barriers") {
        return(.barriers(x, directions, symm, intervalBreaks))
      }
      if(transitionFunction == "areas") {
        return(.areas(x, directions))
      }
    } else {
      if (directions %in% c(4, 8)) {
        if (.isGlobalLonLat(x)) {
          message("Global lat/lon rasters are not supported under new ",
                  "optimizations for 4 and 8 directions with custom transition ",
                  "functions. Falling back to old method.")
          return(.TfromR_old(x, transitionFunction, directions, symm))
        } else {
          return(.TfromR_new(x, transitionFunction, directions, symm))
        }
      } else {
        return(.TfromR_old(x, transitionFunction, directions, symm))
      }
    }
  }
)

# Modified from https://github.com/andrewmarx/samc/blob/1d9973882477180fa90ca7a570c3a0db8cadfbe2/R/internal-functions.R#L332
.tr_vals_simple <- function(data, fun, dir, sym) {

  nrows <- terra::nrow(data)
  ncols <- terra::ncol(data)

  if (sym) {
    result <- numeric(nrows * ncols * dir / 2)
  } else {
    result <- numeric(nrows * ncols * dir)
  }
  index <- 0

  if (dir == 4) {
    dir <- c(2, 4, 6, 8)
  } else if (dir == 8) {
    dir <- c(1:4, 6:9)
  } else if (dir == 16) {
    # TODO 16 directions
    stop("16 directions not reimplemented yet")
  } else {
    stop("Bad directions", call. = FALSE)
  }

  if (sym) dir <- dir[1:(length(dir)/2)]

  for (r in 1:nrows) {
    vals <- terra::focalValues(data, 3, r, 1)

    for (c in 1:ncols) {
      v <- vals[c, 5]

      for (d in dir) {
        index <- index + 1

        if (is.finite(v) & is.finite(vals[c, d])) {
          result[index] <- fun(c(vals[c, d], v))
        }
      }
    }
  }

  result
}

# Modified from https://github.com/andrewmarx/samc/blob/1d9973882477180fa90ca7a570c3a0db8cadfbe2/R/internal-functions.R#L13
.TfromR_new <- function(x, tr_fun, dir, sym) {
  extent <- raster::extent(x)
  crs <- raster::projection(x, asText=FALSE)

  x <- terra::rast(x)

  tr_vals <- .tr_vals_simple(x, tr_fun, dir, sym)

  nedges <- sum(tr_vals != 0)

  nrows <- terra::nrow(x)
  ncols <- terra::ncol(x)
  ncells <- terra::ncell(x)
  cell_nums <- terra::cells(x)

  rm(x)

  if (dir == 4) {
    dir_vec <- c(1:4)
    offsets <- c(-ncols, -1, 1, ncols)
  } else if (dir == 8) {
    dir_vec <- c(1:8)
    offsets <- c(-ncols - 1, -ncols, -ncols + 1, -1, 1, ncols - 1, ncols, ncols + 1)
  } else if (dir == 16) {
    # TODO 16 directions
    stop("16 directions not reimplemented yet")
  } else {
    stop("Bad directions", call. = FALSE)
  }

  dir_sym <- dir

  if (sym) {
    dir_sym <- dir/2
    dir_vec <- dir_vec[1:dir_sym]
    offsets <- offsets[1:dir_sym]
  }


  mat_p <- integer(ncells + 1)
  mat_x <- numeric(nedges)
  mat_i <- integer(nedges)

  index = 0

  for (i in 1:length(tr_vals)) {
    if (tr_vals[i] != 0) {
      index = index + 1

      # TODO modify for rasters that wrap across date line
      d <- ((i - 1) %% dir_sym) + 1
      p2 <- ((i - 1) %/% dir_sym) + 1
      p1 <- p2 + offsets[d]

      mat_p[p2 + 1] <- mat_p[p2 + 1] + 1
      mat_x[index] <- tr_vals[i]
      mat_i[index] <- p1
    }
  }

  for (i in 2:length(mat_p)) {
    mat_p[i] <- mat_p[i] + mat_p[i - 1]
  }

  mat_i <- mat_i - 1

  if (!all(mat_x >= 0)) {
    warning("transition function gives negative values")
  }

  if (sym) {
    mat <- new("dsCMatrix")
  } else {
    mat <- new("dgCMatrix")
  }

  mat@Dim <- c(as.integer(ncells), as.integer(ncells))

  mat@p <- as.integer(mat_p)
  mat@i <- as.integer(mat_i)
  mat@x <- mat_x

  new("TransitionLayer",
      nrows = as.integer(nrows),
      ncols = as.integer(ncols),
      extent = extent,
      crs = crs,
      transitionMatrix = mat,
      transitionCells = 1:ncells,
      matrixValues = "conductance")
}

.TfromR_old <- function(x, transitionFunction, directions, symm) {
  tr <- new("TransitionLayer",
            nrows = as.integer(raster::nrow(x)),
            ncols = as.integer(raster::ncol(x)),
            extent = raster::extent(x),
            crs = raster::projection(x, asText = FALSE),
            transitionMatrix = Matrix(0, raster::ncell(x), raster::ncell(x)),
            transitionCells = 1:raster::ncell(x))

  transitionMatr <- transitionMatrix(tr)
  Cells <- which(!is.na(raster::getValues(x)))
  adj <- raster::adjacent(x, cells = Cells, pairs = TRUE,
                          target = Cells,
                          directions = directions)
  if(symm) { adj <- adj[adj[, 1] < adj[, 2], ] }
  dataVals <- cbind(raster::getValues(x)[adj[, 1]],
                    raster::getValues(x)[adj[, 2]])
  transition.values <- apply(dataVals, 1, transitionFunction)

  if(!all(transition.values >= 0)){
    warning("transition function gives negative values")
  }

  transitionMatr[adj] <- as.vector(transition.values)
  if(symm) {
    transitionMatr <- forceSymmetric(transitionMatr)
  }
  transitionMatrix(tr) <- transitionMatr
  matrixValues(tr) <- "conductance"

  return(tr)
}

.barriers <- function(x, directions, symm, intervalBreaks) {
	Xlayer <- new("TransitionLayer",
		nrows = as.integer(nrow(x)),
		ncols = as.integer(ncol(x)),
		extent = extent(x),
		crs = projection(x, asText = FALSE),
		transitionMatrix = Matrix(0, ncell(x), ncell(x)),
		transitionCells = 1:ncell(x))
	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)

	if(x@data@isfactor) {

		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)

		if(symm) {
			maxn <- (n^2 - n)/2
			for(i in 1:maxn) {
				j <- .matrIndex(i, n)
				XlayerNew <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacent(x, cells = cells1, pairs = TRUE,
				                 target = cells2,
				                 directions = directions)
				adj2 <- adjacent(x, cells = cells2,
				                 pairs = TRUE,
				                 target = cells1,
				                 directions = directions)
				adj <- rbind(adj1, adj2)
				XlayerNew[adj] <- 1
				Xstack <- stack(Xstack, XlayerNew)
			}
		} else {
			maxn <- (n^2 - n)/2
			for(i in 1:maxn) {
				j <- .matrIndex(i,n)
				XlayerNew1 <- Xlayer
				XlayerNew2 <- Xlayer
				cells1 <- which(getValues(x) == vals[j[1]])
				cells2 <- which(getValues(x) == vals[j[2]])
				adj1 <- adjacent(x,
				                 cells = cells1,
				                 pairs = TRUE,
				                 target = cells2,
				                 directions = directions)
				adj2 <- adjacent(x,
				                 cells = cells2,
				                 pairs = TRUE,
				                 target = cells1,
				                 directions = directions)
				XlayerNew1[adj1] <- 1
				XlayerNew2[adj2] <- 1
				Xstack <- stack(Xstack, XlayerNew1, XlayerNew2)
			}
		}
	} else {
		Xmin <- transition(x, min, directions)
		Xmax <- transition(x, max, directions)
		index1 <- adjacent(x,
		                   cells = 1:ncell(x),
		                   pairs = TRUE,
		                   target = 1:ncell(x),
		                   directions = directions)
		XminVals <- Xmin[index1]
		XmaxVals <- Xmax[index1]

		if(symm == TRUE) {
			for(i in 1:length(intervalBreaks)) {
				index2 <- index1[XminVals < intervalBreaks[i] & XmaxVals > intervalBreaks[i], ]
				XlayerNew <- Xlayer
				XlayerNew[index2] <- 1
				Xstack <- stack(Xstack, XlayerNew)
			}
		}

		if(symm=="up" | symm=="down") {
		  stop("not implemented yet")
		}
	}

	Xstack <- Xstack[[2:nlayers(Xstack)]]
	return(Xstack)
}


.areas <- function(x, directions) {

	Xlayer <- new("TransitionLayer",
		nrows = as.integer(nrow(x)),
		ncols = as.integer(ncol(x)),
		extent = extent(x),
		crs = projection(x, asText = FALSE),
		transitionMatrix = Matrix(0, ncell(x), ncell(x)),
		transitionCells = 1:ncell(x))

	matrixValues(Xlayer) <- "resistance"
	Xstack <- as(Xlayer, "TransitionStack") * 0
	#Xstack@transition <- vector(list,...)

	if(x@data@isfactor) {
		vals <- unlist(x@data@attributes[[1]])
		n <- length(vals)

		for(i in 1:n) {
			transitionFunction <- function(v) { return(sum(v == i) / 2) }
			XlayerNew <- .TfromR_old(x, transitionFunction, directions, symm = TRUE) # TODO integrate with new optimizations
			Xstack <- stack(Xstack, XlayerNew)
		}
	} else {
		warning("not yet implemented for raster with non-factor",
		        " variables. Contact author.")
	}

	Xstack <- Xstack[[2:nlayers(Xstack)]]

	return(Xstack)
}


#' @export
setMethod(
  "transition",
  signature(x = "RasterBrick"),
  def = function(x, transitionFunction = "mahal", directions) {
		if (transitionFunction != "mahal") {
			stop("only Mahalanobis distance method implemented for RasterBrick \n")
		}

		xy <- cbind(1:ncell(x), getValues(x))
		xy <- na.omit(xy)

		dataCells <- xy[, 1]

		adj <- adjacent(x, cells = dataCells, pairs = TRUE,
		                target = dataCells, directions = directions)

		x.minus.y <- raster::getValues(x)[adj[, 1], ] - raster::getValues(x)[adj[, 2], ]

		cov.inv <- solve(cov(xy[, -1]))

		mahaldistance <- apply(x.minus.y, 1, function(x) { sqrt((x%*%cov.inv)%*%x) })
		mahaldistance <- mean(mahaldistance)/(mahaldistance + mean(mahaldistance))

		transitiondsC <- new("dsCMatrix",
				p = as.integer(rep(0, ncell(x) + 1)),
				Dim = as.integer(c(ncell(x), ncell(x))),
				Dimnames = list(as.character(1:ncell(x)), as.character(1:ncell(x))))

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
