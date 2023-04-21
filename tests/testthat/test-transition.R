context("transition")

# TODO implement lonlat tests
tr_params = data.frame(
  sym = c(TRUE, FALSE)
)

for (test in testlist) {
  for (i in 1:nrow(tr_params)) {
    tr = transition(test, function(x) 1/mean(x), 8, symm = tr_params$sym[i])
    trmat = transitionMatrix(tr)
    geotr = geoCorrection(tr)
    geotrmat = transitionMatrix(geotr)
    
    cells = which(!is.na(raster::getValues(test)))
    adj = raster::adjacent(test, cells = cells, pairs = TRUE, target = cells, directions = 8)
    tr_count = nrow(adj)
    if (tr_params$sym[i]) tr_count = tr_count/2
    

    # Validate TransitionLayer attributes
    test_that("TransitionLayer attributes", {
      expect_equal(nrow, nrow(tr))
      expect_equal(ncol, ncol(tr))
      expect_equal(ncell, ncell(tr))
      expect_equal(ext, extent(tr))
    })
    
    
    # Validate TransitionLayer matrix
    test_that("TransitionLayer matrix", {
      expect_equal(dim(trmat), c(ncell, ncell))
      expect_equal(length(trmat@i), tr_count)
      expect_equal(length(trmat@x), tr_count)
      expect_equal(trmat@p[length(trmat@p)], tr_count)
      # TODO actual values
    })
    
    
    # Validate geocorrected TransitionLayer matrix
    test_that("TransitionLayer geocorrected matrix", {
      expect_equal(dim(geotrmat), c(ncell, ncell))
      expect_equal(length(geotrmat@i), tr_count)
      expect_equal(length(geotrmat@x), tr_count)
      expect_equal(geotrmat@p[length(geotrmat@p)], tr_count)
      # TODO actual values
    })
  }
}
