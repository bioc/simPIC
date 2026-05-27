# Assuming you have already loaded counts from the RDS file
counts <- readRDS(system.file("extdata", "test.rds", package = "simPIC"))

test_that("simPICestimate works with dgCMatrix input", {
  library(Matrix)
  
  # Ensure counts is a valid dgCMatrix object (if it's not already)
  mat <- as(counts, "dgCMatrix")
  
  simPIC_obj <- newsimPICcount()
  
  # Run estimate
  est_obj <- simPICestimate(mat, object = simPIC_obj, pm.distr = "weibull", method = "single", verbose = FALSE)
  
  # Check output class
  expect_s4_class(est_obj, "simPICcount")
  
  # Check that some expected slots are populated
  expect_true(length(est_obj@nPeaks) > 0)
})

test_that("simPICestimate works with SingleCellExperiment input", {
  library(SingleCellExperiment)
  
  # Simulate SCE using counts
  sce <- SingleCellExperiment(assays = list(counts = counts))
  
  est_obj <- simPICestimate(sce, pm.distr = "weibull", method = "single", verbose = FALSE)
  
  expect_s4_class(est_obj, "simPICcount")
  expect_true(length(est_obj@nPeaks) > 0)
})

test_that("simPICestimate works with weibull distribution", {
  object <- simPICestimate(counts, pm.distr = "weibull")
  expect_true(validObject(object))
})

test_that("simPICestimate works with gamma distribution", {
  object <- simPICestimate(counts, pm.distr = "gamma")
  expect_true(validObject(object))
})

test_that("simPICestimate estimates pareto parameters with verbose FALSE", {
  object <- simPICestimate(counts, pm.distr = "pareto", verbose = FALSE)
  defaults <- newsimPICcount()

  expect_true(validObject(object))
  expect_identical(object@pm.distr, "pareto")
  expect_false(isTRUE(all.equal(object@peak.mean.shape, defaults@peak.mean.shape)))
  expect_false(isTRUE(all.equal(object@peak.mean.scale, defaults@peak.mean.scale)))
})

test_that("simPICestimate estimates lngamma parameters with verbose FALSE", {
  object <- simPICestimate(counts, pm.distr = "lngamma", verbose = FALSE)
  defaults <- newsimPICcount()

  expect_true(validObject(object))
  expect_identical(object@pm.distr, "lngamma")
  expect_false(isTRUE(all.equal(object@peak.mean.pi, defaults@peak.mean.pi)))
  expect_false(isTRUE(all.equal(object@peak.mean.meanlog, defaults@peak.mean.meanlog)))
  expect_false(isTRUE(all.equal(object@peak.mean.sdlog, defaults@peak.mean.sdlog)))
})

test_that("simPICestimate works with SingleCellExperiment", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  object <- simPICestimate(sce, pm.distr = "weibull")
  expect_true(validObject(object))
})

test_that("simPICestimate SingleCellExperiment honours verbose FALSE", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  object <- simPICestimate(sce, pm.distr = "pareto", verbose = FALSE)

  expect_identical(object@pm.distr, "pareto")
  expect_true(validObject(object))
})

test_that("simPICestimate works with SingleCellExperiment without counts", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(TEST = counts)
  )
  expect_warning(simPICestimate(sce, pm.distr = "weibull"), "counts assay is missing")
})
