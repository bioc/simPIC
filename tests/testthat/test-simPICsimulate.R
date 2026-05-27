# Checks simPICsimulate

library(testthat)

test_that("simPICsimulatesingle returns a SingleCellExperiment", {
  skip_if_not_installed("SingleCellExperiment")
  
  sim <- simPICsimulatesingle()
  
  # Check class
  expect_s4_class(sim, "SingleCellExperiment")
  
  # Check dimensions (some defaults should be set in newsimPICcount())
  expect_true(ncol(sim) > 0)
  expect_true(nrow(sim) > 0)
  
  # Check colData and rowData are not empty
  expect_true(ncol(colData(sim)) >= 1)
  expect_true(ncol(rowData(sim)) >= 1)
  
  # Check counts assay exists
  expect_true("counts" %in% SummarizedExperiment::assayNames(sim))
})

test_that("simPICsimulatemulti assigns groups when nGroups > 1", {
  custom_obj <- newsimPICcount(nGroups = 3, group.prob=c(0.5,0.25,0.25))  # force multiple groups
  sim <- simPICsimulatemulti(object = custom_obj)
  
  expect_s4_class(sim, "SingleCellExperiment")
  expect_true("Group" %in% colnames(colData(sim)))
  expect_gt(length(unique(colData(sim)$Group)), 1)
})


test_that("simPICsimulate function works correctly with weibull", {
  sim_result <- simPICsimulate()
  expect_true(inherits(sim_result, "SingleCellExperiment"))
})

test_that("simPICsimulate function works correctly with gamma", {
  sim_result <- simPICsimulate(pm.distr = "gamma")
  expect_true(inherits(sim_result, "SingleCellExperiment"))
})

test_that("simPICsimulate function works correctly with lngamma", {
  params <- newsimPICcount()
  sim_result <- simPICsimulate(pm.distr = "lngamma")
  expect_true(inherits(sim_result, "SingleCellExperiment"))
})

test_that("simPICsimulate function works correctly with pareto", {
  params <- newsimPICcount()
  sim_result <- simPICsimulate(pm.distr = "pareto")
  expect_true(inherits(sim_result, "SingleCellExperiment"))
})

test_that("simPICsimulate returns counts as the first assay", {
  sim_result <- simPICsimulate(verbose = FALSE)
  expect_identical(SummarizedExperiment::assayNames(sim_result)[1], "counts")
})

test_that("simPICsimulate stores the selected peak-mean distribution", {
  sim_result <- simPICsimulate(pm.distr = "lngamma", verbose = FALSE)
  expect_identical(S4Vectors::metadata(sim_result)$Params@pm.distr, "lngamma")
})

test_that("nGroups is 1, switching to default mode", {
  params <- newsimPICcount()
  expect_warning(
    simPICsimulate(params,
                  method = "groups",
                  group.prob = c(1)
    ),
    "nGroups is 1, switching to default mode"
  ) 
  
})

test_that("simPICsimulate throws error on invalid input", {
  expect_error(simPICsimulate(object = "not_a_simPICcount"))
})
