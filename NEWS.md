## Version 1.5.3 (2025-04-28)
* Extended `simPICcountclass` object to include: nGroups, batch, differential
  accessibility, and bcv parameters.
* Updated `estimateBCV` function in `simPICestimate`.
* Updates in `simPICsimulate`:
  - `simPICsimBatchEffects`
  - `simPICsimBatchCellMeans`
  - `simPICsimulatemultiDA`
  - `simPICsimulateBCVmeans`
  - `simPICsimulateTrueCountsGroups`
* Updated vignette to demonstrate new functionality--- multiple cell-types
* Removed redundancy in citation in `DESCRIPTION`

## Version 1.5.2 (2025-04-21)
* Added package logo to vignette.

## Version 1.5.1 (2025-04-21)
* Fixing typo in vignette, changed `lognormal` to `lognormal-gamma`

## Version 0.99.7 (2024-04-14)
* Updating title.
* Simplifying estimate sparsity description.

## Version 0.99.6 (2024-03-16)
* Replacing sapply with vapply in `simPICsimulate` line 238.

## Version 0.99.5 (2024-03-16)
* Addressing notes after reviewer comments.
* Updated simPIC count documentation. 
* Updated 1:nCells to seq_len in `simPICsimulate` line 239.
* Added documentation for testdata.R.

## Version 0.99.4 (2024-03-10)
* Fixed issue with system files found that should not be Git
  tracked in BiocCheck::BiocCheckGitClone

## Version 0.99.3 (2024-03-10)
* Updating R dependency to R (>=4.4.0)

## Version 0.99.2 (2024-03-09)
* Fixing issue with system files found in BiocCheck

## Version 0.99.1 (2024-03-09)
* Addressed reviewer suggestions
* Major changes include improved docs, including a package man page, text file 
and code to reproduce test data

## Version 0.99.0 (2024-02-28)

* Submitted to Bioconductor