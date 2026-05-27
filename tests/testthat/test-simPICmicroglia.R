test_that("packaged Microglia example files load", {
  skip_if_not_installed("VariantAnnotation")

  paths <- simPIC:::simPICmicrogliaExamplePaths()

  expect_true(all(file.exists(unlist(paths))))

  sce <- readRDS(paths$sce)
  expect_s4_class(sce, "SingleCellExperiment")
  expect_true(all(c("Sample", "Library", "Batch", "sample_batch", "Pathology") %in%
    colnames(SummarizedExperiment::colData(sce))))

  sample.map <- utils::read.delim(paths$sample_map, stringsAsFactors = FALSE)
  expect_true(all(c("vcf_id", "Sample") %in% colnames(sample.map)))

  peak.annot <- utils::read.delim(paths$peak_annot, stringsAsFactors = FALSE)
  expect_true(all(c("seqnames", "start", "end", "peak_id") %in%
    colnames(peak.annot)))
})

test_that("simPICMicrogliaExample uses bigcounts and aggregated means", {
  skip_if_not_installed("splatter")
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("scater")

  out <- suppressWarnings(simPICMicrogliaExample(
    min.cells = 20,
    eqtl.n = 0.05,
    similarity.scale = 1.5,
    eqtl.dist = 1e8,
    pca.ntop = 100,
    pca.components = 3,
    sparsify = FALSE,
    seed = 1,
    verbose = FALSE
  ))

  expect_s4_class(out$sce, "SingleCellExperiment")
  expect_s4_class(out$bigcounts, "SingleCellExperiment")
  expect_s4_class(out$aggregated, "SingleCellExperiment")
  expect_s4_class(out$params, "SplatPopParams")
  expect_s4_class(out$sim, "SingleCellExperiment")
  expect_false(splatter::getParam(out$params, "pop.quant.norm"))
  expect_equal(
    ncol(out$aggregated),
    length(unique(out$sce$sample_batch))
  )
  expect_equal(
    ncol(out$bigcounts),
    max(as.integer(table(out$sce$sample_batch)))
  )
  expect_lte(length(out$plot_samples), 4)
  expect_gte(length(out$plot_samples), 1)
  expect_true(all(c("real", "simulated") %in% names(out$plots)))
  expect_true("comparison" %in% names(out$plots))
  expect_s3_class(out$plots$real$cell_pca_by_sample, "ggplot")
  expect_s3_class(out$plots$simulated$cell_pca_by_sample, "ggplot")
  expect_s3_class(out$plots$comparison$cell_pca_by_sample, "ggplot")
})
