# Checks simPICestimate

counts <- readRDS(system.file("extdata", "test.rds",
  package = "simPIC"
))

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
