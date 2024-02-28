    # Checks simPICsimulate
    # library(testthat)
    #
    # test.object <- newsimPICcount(nPeaks = 10000, nCells = 200,
    #                               lib.size.meanlog = 9.04, lib.size.sdlog = 0.7916,
    #                               peak.mean.shape = 0.654, peak.mean.scale = 0.1194,
    #                               sparsity = runif(200, min = 0.9, max = 0.99))
    #
    # test_that("simPIC simulate works", {
    #   expect_true(validObject(simPICsimulate(test.object, pm.distr = "weibull", seed =34167)))
    #
    # })


    library(testthat)

    # Define a test for simPICsimulate function
    test_that("simPICsimulate function works correctly with weibull", {
      # Define parameters for simPICcount object
      params <- newsimPICcount()

      # Create a simPICcount object
      # simPIC_object <- params

      # Run simPICsimulate function
      sim_result <- simPICsimulate(
        object = params, verbose = FALSE,
        pm.distr = "weibull"
      )

      # Check if the result is a SingleCellExperiment object
      expect_true(inherits(sim_result, "SingleCellExperiment"))
    })
    #   # Check if the dimensions of the result match the expected values
    #   expect_equal(assays(sim_result)$counts@Dim)
    # })

    # Define a test for simPICsimulate function
    test_that("simPICsimulate function works correctly with gamma", {
      # Define parameters for simPICcount object
      params <- newsimPICcount()

      # Create a simPICcount object
      # simPIC_object <- params

      # Run simPICsimulate function
      sim_result <- simPICsimulate(
        object = params, verbose = FALSE,
        pm.distr = "gamma"
      )

      # Check if the result is a SingleCellExperiment object
      expect_true(inherits(sim_result, "SingleCellExperiment"))
    })

    # Define a test for simPICsimulate function
    test_that("simPICsimulate function works correctly with gamma", {
      # Define parameters for simPICcount object
      params <- newsimPICcount()

      # Create a simPICcount object
      # simPIC_object <- params

      # Run simPICsimulate function
      sim_result <- simPICsimulate(
        object = params, verbose = FALSE,
        pm.distr = "lngamma"
      )

      # Check if the result is a SingleCellExperiment object
      expect_true(inherits(sim_result, "SingleCellExperiment"))
    })
