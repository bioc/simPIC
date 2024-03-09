# Checks simPICsimulate

library(testthat)

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
