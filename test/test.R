source("./R/fastbmdR_main.R")
source("./R/fastbmdR_utils.R")
source("./R/example_data.R")

library(drc)
library(parallel)
library(data.table)

set.seed(42)

test_res <- readRDS("./test/test.rds")

data(example_data) # this creates a dataframe called example_data

models <- c("Exp2", "Exp3", "Exp4", "Exp5", "Poly2", "Lin", "Power", "Hill")
dose <- c(rep(0, 5), rep(25, 5), rep(100, 5),
          rep(200, 5), rep(300, 5), rep(400, 5))
ncpus <- 1

fit_obj <- PerformCurveFitting(data = example_data,
                               dose = dose,
                               ncpus = ncpus,
                               models = models)

fit_obj <- FilterDRFit(fit_obj, lof.pval = 0.1, filt.var = "AIC.model")

bmd_res <- PerformBMDCalc(fit_obj,
                          ncpus = ncpus,
                          num.sds = 1,
                          bmr.method = "sample.mean")

bmd_pass <- bmd_res[bmd_res$all.pass, ]

print("Does the package still produce the same output?")
print(all.equal(bmd_pass, test_res))