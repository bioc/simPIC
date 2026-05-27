test_that("splatPopEstimatePeak estimates params from SingleCellExperiment", {
  skip_if_not_installed("splatter")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(101)
  gene_means <- rgamma(60, shape = 2, rate = 0.4)
  cell_scales <- runif(48, min = 0.7, max = 1.4)
  counts <- vapply(
    cell_scales,
    function(scale) {
      rnbinom(60, mu = gene_means * scale, size = 0.2)
    },
    numeric(60)
  )
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = data.frame(
      Sample = rep(paste0("S", 1:12), each = 4),
      Batch = rep(c("B1", "B2"), each = 24)
    )
  )

  params <- splatPopEstimatePeak(
    counts = sce,
    sample.col = "Sample",
    batch.col = "Batch",
    aggregate.by = "sample",
    min.cells = 2,
    pop.cv.bins = 10,
    verbose = FALSE
  )

  expect_s4_class(params, "SplatPopParams")

  expected <- edgeR::estimateDisp(counts, design = matrix(1, ncol(counts), 1))
  expect_equal(
    params@bcv.common,
    -0.3 + 0.15 * expected$common.dispersion,
    tolerance = 1e-8
  )
  expect_equal(params@bcv.df, expected$prior.df, tolerance = 1e-8)
})

test_that("splatPopEstimatePeak requires explicit means for matrix input", {
  skip_if_not_installed("splatter")

  counts <- matrix(rpois(20 * 12, lambda = 3), nrow = 20, ncol = 12)

  expect_error(
    splatPopEstimatePeak(counts = counts, verbose = FALSE),
    "Automatic mean aggregation requires"
  )
})

test_that("splatPopSimulatePeak simulates counts from peak annotations", {
  skip_if_not_installed("splatter")
  skip_if_not_installed("VariantAnnotation")
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("IRanges")

  params <- splatter::newSplatPopParams(nGenes = 20)
  params <- splatter::setParams(params, batchCells = c(20))

  gff <- splatter::mockGFF()
  gff <- gff[gff[[3]] == "gene", ][seq_len(20), ]
  starts <- as.integer(as.character(gff[[4]]))
  ends <- as.integer(as.character(gff[[5]]))
  peak.gr <- GenomicRanges::GRanges(
    seqnames = gff[[1]],
    ranges = IRanges::IRanges(
      start = pmin(starts, ends),
      end = pmax(starts, ends)
    )
  )
  names(peak.gr) <- paste0("peak_", seq_len(length(peak.gr)))

  sim <- splatPopSimulatePeak(
    params = params,
    vcf = splatter::mockVCF(),
    gff = peak.gr,
    sparsify = FALSE,
    verbose = FALSE
  )

  expect_s4_class(sim, "SingleCellExperiment")
  expect_true("counts" %in% SummarizedExperiment::assayNames(sim))
  expect_identical(
    S4Vectors::metadata(sim)$simPIC.population.model,
    "splatPopPeak"
  )
})
