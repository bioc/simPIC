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

test_that("simPICestimate works with SingleCellExperiment", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )
  object <- simPICestimate(sce, pm.distr = "weibull")
  expect_true(validObject(object))
})

test_that("simPICestimate works with SingleCellExperiment without counts", {
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(TEST = counts)
  )
  expect_warning(simPICestimate(sce, pm.distr = "weibull"), "counts assay is missing")
})
