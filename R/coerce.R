#' @importClassesFrom raster RasterLayer
#' @importClassesFrom raster RasterBrick
#' @importClassesFrom Matrix sparseMatrix
#' @exportMethod coerce
setAs("TransitionLayer", "sparseMatrix", function(from) {
  from@transitionMatrix
})

setAs("TransitionData", "sparseMatrix", function(from) {
  from@transitionMatrix
})

setAs("TransitionLayer", "RasterLayer", function(from) {
  raster(
    xmn = xmin(from),
    xmx = xmax(from),
    ymn = ymin(from),
    ymx = ymax(from),
    nrows = nrow(from),
    ncols = ncol(from),
    crs = projection(from)
  )
})

setAs("RasterLayer", "TransitionLayer", function(from) {
  new(
    "TransitionLayer",
    nrows = as.integer(from@nrows),
    ncols = as.integer(from@ncols),
    extent = extent(c(
        xmin = from@xmin,
        xmax = from@xmax,
        ymin = from@ymin,
        ymax = from@ymax
    )),
    crs = raster::projection(from, asText = FALSE)
  )
})

setAs("TransitionLayer", "TransitionStack", function(from) {
  TS <- methods::new(
    "TransitionStack",
    nrows = as.integer(nrow(from)),
    ncols = as.integer(ncol(from)),
    extent = raster::extent(c(
      xmin = xmin(from),
      xmax = xmax(from),
      ymin = ymin(from),
      ymax = ymax(from)
    )),
    crs = raster::projection(from, asText = FALSE)
  )
  TS@transition <- list(methods::as(from, "TransitionData"))
  return(TS)
})

setAs("TransitionLayer", "TransitionData", function(from) {
  TD <- new(
    "TransitionData",
    transitionMatrix = from@transitionMatrix,
    transitionCells = from@transitionCells,
    matrixValues = from@matrixValues
  )
  return(TD)
})
