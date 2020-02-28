## ----setup, include=FALSE-----------------------------------------------------
TRAVIS <- !identical(tolower(Sys.getenv("TRAVIS")), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  purl = TRAVIS
)

## ----fig0, echo=FALSE, message=FALSE, fig.cap="Fig.1 : Three commonly used ways to generate neighbour graphs for distance calculations.", fig.align="center"----
library("gdistance")
rex <- raster(matrix(1,4,4))

a <- rep(c(1.3333), times=5)
b <- c(-1.3333, -0.6666, 0, 0.6666, 1.3333)

x1 <- c(-a, b)
x2 <- c(a, b)
y1 <- c(b, -a)
y2 <- c(b, a)
x <- cbind(x1,x2)
y <- cbind(y1,y2)

par(mfrow=c(1,3), mar= c(2,2,2,2), oma = c(0,0,0,0) + 0.1, cex.main=1)

x4 <- transition(rex, mean, 4)
g4 <- graph.adjacency(transitionMatrix(x4), mode="undirected")
gridLayout <- xyFromCell(x4, 1:ncell(x4))
plot(g4,layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="4 neighbours")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g4, layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

x8 <- transition(rex, mean, 8)
g8 <- graph.adjacency(transitionMatrix(x8), mode="undirected")
plot(g8,layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="8 neighbours")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g8, layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

x16 <- transition(rex, mean, 16)
g16 <- graph.adjacency(transitionMatrix(x16), mode="undirected")
plot(g16, layout=gridLayout, edge.color="black", vertex.color="black", vertex.label=NA, main="16 neighbours")
for(i in 1:dim(x)[1]){lines(x[i,],y[i,], col="lightgray")}
plot(g16,layout=gridLayout, add=TRUE, edge.color="black", vertex.color="black", vertex.label=NA)

## ----gdistance-1--------------------------------------------------------------
library("gdistance")
set.seed(123)
r <- raster(ncol=3,nrow=3)
r[] <- 1:ncell(r)
r

## ----figure1, fig=TRUE, height = 3.9, echo=FALSE, fig.cap="Fig. 2: Lon-lat grid with cell numbers of a 3 × 3 raster", fig.align="center"----
plot(r, main="r", xlab="Longitude (degrees)", ylab="Latitude (degrees)")
text(r)

## ----gdistance-3--------------------------------------------------------------
r[] <- 1
tr1 <- transition(r, transitionFunction=mean, directions=8)

## ----gdistance-4--------------------------------------------------------------
tr1

## ----gdistance-5--------------------------------------------------------------
r[] <- runif(9)
ncf <- function(x) max(x) - x[1] + x[2]
tr2 <- transition(r, ncf, 4, symm=FALSE)
tr2

## ----gdistance-6--------------------------------------------------------------
tr3 <- tr1*tr2
tr3 <- tr1+tr2
tr3 <- tr1*3
tr3 <- sqrt(tr1)

