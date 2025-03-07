# fastbmdR

`fastbmdR` is an R package designed for performing computationally-efficient dose-response analysis of omics data. It's intial implementation in our [2019 publication in Bioinformatics](https://academic.oup.com/bioinformatics/article/37/7/1035/5881631). The fastbmd curve fitting functions are hosted in the [ExpressAnalyst webtool](https://www.expressanalyst.ca), and you can follow Basic Protocol 10 in [this publication](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/cpz1.922) to learn how to use ExpressAnalyst for dose-response analysis.

## Installation

You can install the `fastbmdR` package directly from GitHub using the `devtools` package. To install `devtools`, run the following command in your R console:

```r
install.packages("devtools")
```

Once devtools is installed, you can install fastbmdR using:
```r
devtools::install_github("jessica-ewald/fastbmdR")
```

After installation, load the package using:
```r
library(fastbmdR)
```

## Usage

Load the example data and specify basic parameters:
```r
data("example_data")
models <- c("Exp2", "Exp3", "Exp4", "Exp5", "Poly2", "Lin", "Power", "Hill")
dose <- c(rep(0, 5), rep(25, 5), rep(100, 5), rep(200, 5), rep(300, 5), rep(400, 5))
ncpus <- 1
```
First, perform the initial curve fitting.
```r
fitObj <- PerformCurveFitting(data = example_data, dose = dose, ncpus = ncpus, models = models)
```
Next, filter to select the best curve fit for each gene. All curve fits must have a lack-of-fit p-value greater than 0.1 and we are using the model AIC values to select the best fit.
```r
fitObj <- FilterDRFit(fitObj, lof.pval = 0.1, filt.var = "AIC.model")
```
Finally, calculate the benchmark dose (BMD) for each curve fit. We also want to filter out low quality BMDs by only keeping those with the "all.pass" flag.
```r
bmd.res <- PerformBMDCalc(fitObj, ncpus = ncpus, num.sds = 1, bmr.method = "sample.mean")
bmd.pass <- bmd.res[bmd.res$all.pass, ]
```

## Testing
If changes are made to the fastbmdR_main.R or fastbmdR_utils.R, run `Rscript test/test.R` from the root fastbmdR directory.
