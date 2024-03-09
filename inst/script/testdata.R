# Script to reproduce default parameters

options(max.print = 30)

counts <- readRDS("inst/extdata/test.rds")
new <- newsimPICcount()
new <- simPICestimate(counts,pm.distr ="weibull")
print(new)
