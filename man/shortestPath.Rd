% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shortestPath.R
\name{shortestPath}
\alias{shortestPath}
\alias{shortestPath,TransitionLayer,Coords,Coords-method}
\title{Shortest path}
\usage{
shortestPath(x, origin, goal, ...)
}
\arguments{
\item{x}{\code{TransitionLayer} object}

\item{origin}{SpatialPoints, vector or matrix with coordinates,
at the moment only the first cell is taken into account}

\item{goal}{SpatialPoints, vector or matrix with coordinates}

\item{...}{Additional argument: output}
}
\value{
Transition object
}
\description{
Calculates the shortest path from an origin to a goal
}
\details{
As an additional argument output, the desired output object can be specified:
 either \dQuote{TransitionLayer} (default),
 \dQuote{TransitionStack} or \dQuote{SpatialLines}.

If there is more than one path either (1) transition values in the
 TransitionLayer get values of 1 / number of paths or (2) the SpatialLines
 object will contain more than one line.
}
\examples{
#example equivalent to that in the documentation on r.cost/r.drain in GRASS
r <- raster(nrows=6, ncols=7, xmn=0, xmx=7, ymn=0, ymx=6, crs="+proj=utm +units=m")

r[] <- c(2, 2, 1, 1, 5, 5, 5,
         2, 2, 8, 8, 5, 2, 1,
         7, 1, 1, 8, 2, 2, 2,
         8, 7, 8, 8, 8, 8, 5,
         8, 8, 1, 1, 5, 3, 9,
         8, 1, 1, 2, 5, 3, 9)

tr <- transition(r, function(x) 1/mean(x), 8)
# 1/mean: reciprocal to get permeability
tr <- geoCorrection(tr)

c1 <- c(5.5,1.5)
c2 <- c(1.5,5.5)

#make a SpatialLines object for visualization
sPath1 <- shortestPath(tr, c1, c2, output="SpatialLines")
plot(r)
lines(sPath1)

#make a TransitionLayer for further calculations
sPath2 <- shortestPath(tr, c1, c2)

plot(raster(sPath2))

}
\seealso{
\code{\link{costDistance}}, \code{\link{accCost}}
}
\author{
Jacob van Etten
}
