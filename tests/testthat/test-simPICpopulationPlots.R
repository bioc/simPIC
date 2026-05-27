test_that("simPICplotPopulationPCA builds PCA plots from matrix input", {
  skip_if_not_installed("scater")

  set.seed(11)
  counts <- matrix(rpois(80 * 24, lambda = 3), nrow = 80, ncol = 24)
  colnames(counts) <- paste0(
    rep(paste0("Sample", 1:6), each = 4),
    "#Cell",
    seq_len(ncol(counts))
  )
  rownames(counts) <- paste0("Peak", seq_len(nrow(counts)))

  out <- suppressWarnings(simPICplotPopulationPCA(
    counts,
    pca.ntop = 50,
    pca.components = 5,
    verbose = FALSE
  ))

  expect_s4_class(out$cell_sce, "SingleCellExperiment")
  expect_s4_class(out$sample_sce, "SingleCellExperiment")
  expect_true(all(c("cell_pca_by_sample", "cells_per_sample", "sample_pca") %in%
    names(out$plots)))
  expect_s3_class(out$plots$cell_pca_by_sample, "ggplot")
  expect_equal(length(unique(out$cell_sce$Sample)), 6)
})

test_that("simPICplotBlusterComparison builds Poptrial-style panels", {
  skip_if_not_installed("scater")
  skip_if_not_installed("bluster")

  set.seed(21)
  counts.real <- matrix(rpois(60 * 18, lambda = 4), nrow = 60, ncol = 18)
  counts.sim <- matrix(rpois(60 * 18, lambda = 4), nrow = 60, ncol = 18)
  rownames(counts.real) <- rownames(counts.sim) <- paste0("Peak", seq_len(60))

  real <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts.real),
    colData = S4Vectors::DataFrame(
      Sample = rep(paste0("Sample", 1:3), each = 6),
      Batch = rep("Library5", 18)
    )
  )
  sim <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts.sim),
    colData = S4Vectors::DataFrame(
      Sample = rep(paste0("Sample", 1:3), each = 6)
    )
  )

  out <- suppressWarnings(simPICplotBlusterComparison(
    real.sce = real,
    simulated.sce = sim,
    sample.col = "Sample",
    batch.col = "Batch",
    subset.batch = "Library5",
    pca.ntop = 40,
    pca.components = 3,
    verbose = FALSE
  ))

  expect_s4_class(out$real_sce, "SingleCellExperiment")
  expect_s4_class(out$simulated_sce, "SingleCellExperiment")
  expect_equal(sort(out$plot_samples), paste0("Sample", 1:3))
  expect_s3_class(out$plots$cell_pca_by_sample, "ggplot")
  expect_s3_class(out$plots$silhouette_width, "ggplot")
  expect_s3_class(out$plots$neighborhood_purity, "ggplot")
})
