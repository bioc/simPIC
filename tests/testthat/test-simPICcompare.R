# check compare function

sim1 <- simPICsimulate(batchCells = c(700))
sim2 <- simPICsimulate(nPeaks = 10000, nCells = 2000, batchCells = c(2000))

test_that("simPICcompare works", {
  skip_if_not_installed("ggplot2")

  comparison <- simPICcompare(list(real = sim1, simPIC = sim2))
  expect_length(comparison, 4)
  expect_true(all(c("RowData", "ColData", "KS", "Plots") %in%
    names(comparison)))
  checkmate::expect_class(comparison$ColData, "data.frame")
  checkmate::expect_class(comparison$RowData, "data.frame")
  checkmate::expect_class(comparison$KS, "data.frame")
  expect_true(all(c(
    "Metric", "Reference", "Dataset", "Statistic", "PValue"
  ) %in% colnames(comparison$KS)))
  expect_equal(sort(unique(comparison$KS$Metric)),
    c("CellSparsity", "LibrarySizes", "PeakMeans"))
  expect_length(comparison$Plots, 8)
  expect_true(all(c(
    "Means", "Variances", "MeanVar", "LibrarySizes",
    "ZerosPeak", "ZerosCell", "MeanZeros", "PeakMeanNzp"
  ) %in%
    names(comparison$Plots)))
  for (plot in names(comparison$Plots)) {
    checkmate::expect_class(comparison$Plots[[plot]], "ggplot")
  }
})
