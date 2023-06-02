context("transition")

fun_list = list(
  function(x) 1/mean(x), # Symmetric
  function(x) max(x) - x[1] + x[2] # Asymmetric
)

# TODO implement lonlat tests
tr_params = expand.grid(
  directions = c(4, 8, 16),
  sym = c(TRUE, FALSE),
  fun = 1:length(fun_list)
)

for (test in testlist) {
  for (i in 1:nrow(tr_params)) {
    tr = transition(test, fun_list[[tr_params$fun[i]]], directions = tr_params$directions[i], symm = tr_params$sym[i])
    tr_old = .TfromR_old(test, fun_list[[tr_params$fun[i]]], directions = tr_params$directions[i], symm = tr_params$sym[i])

    trmat = transitionMatrix(tr)
    geotr = geoCorrection(tr)
    geotrmat = transitionMatrix(geotr)

    cells = which(!is.na(raster::getValues(test)))
    adj = raster::adjacent(test, cells = cells, pairs = TRUE, target = cells, directions = tr_params$directions[i])
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

      expect_identical(trmat, transitionMatrix(tr_old))

      # TODO compare against manually calculated values
    })


    # Validate geocorrected TransitionLayer matrix
    test_that("TransitionLayer geocorrected matrix", {
      expect_equal(dim(geotrmat), c(ncell, ncell))
      expect_equal(length(geotrmat@i), tr_count)
      expect_equal(length(geotrmat@x), tr_count)
      expect_equal(geotrmat@p[length(geotrmat@p)], tr_count)

      # TODO compare against manually calculated values
    })
  }
}
