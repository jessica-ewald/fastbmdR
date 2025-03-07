##### Other functions for fastbmdR ######

### Hill model and starting values
formHill <- as.formula(signal ~ c + (d - c) / (1 + (dose / e)^b))
startvalHillnls2 <- function(x, y, xm, ym, increase)
# requires the definition of increase from min and max values
# which is the first one
# inputs
# - x values of the dose
# - y values the corresponding signal
# - xm unique values of the dose (sorted by dose)
# - ym means of the signal at each value of xm (sorted by dose)
{
  maxi <- max(y, na.rm = TRUE)
  mini <- min(y, na.rm = TRUE)
  ampl <- maxi - mini

  # inflate maxi and mini so all values are inside [mini; maxi]
  maxi <- maxi + 0.001 * ampl
  mini <- mini - 0.001 * ampl

  # initial value of c
  c <- ifelse(increase, maxi, mini)
  # initial value of d
  d <- ifelse(increase, mini, maxi)
  # initial value of e and b from regression
  yreg <- log((d - c) / (y[x != 0] - c) - 1)
  xreg <- log(x[x != 0])
  reg <- lm(yreg ~ xreg)
  b <- reg$coefficients[2]
  e <- reg$coefficients[1] / (-b)
  startval <- list(b = b, c = c, d = d, e = e)
}


### Exp2 model
# define the model formula
formExp2 <- as.formula(signal ~ e * exp(b * dose))
# get starting values
startvalExp2 <- function(xm, ym)
# inputs
# - xm unique values of the dose (sorted by dose)
# - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of a
  e <- ym[1]

  # transform y for regression
  yreg <- log(ym[xm != 0])
  reg <- lm(yreg ~ xm[xm != 0])

  # estimate slope from regression
  b <- coef(reg)[2]

  startval <- list(e = e, b = b)
}


### Exp3 model
# define the model formula
formExp3 <- as.formula(signal ~ e * (exp(sign(b) * (abs(b) * dose)^d)))

# get starting values
startvalExp3 <- function(xm, ym)
# inputs
# - xm unique values of the dose (sorted by dose)
# - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of e
  e <- ym[1]

  # transform y for regression
  yreg <- log(ym[xm != 0])
  reg <- lm(yreg ~ xm[xm != 0])

  # estimate b and d from regression
  b <- coef(reg)[2]

  d <- (exp(coef(reg)[1])) / e

  startval <- list(e = e, b = b, d = d)
}


### Exp4 model
formExp4 <- as.formula(signal ~ e * (c - (c - 1) * exp((-1) * b * dose)))

# get starting values
startvalExp4 <- function(xm, ym, ad.dir)
# inputs
# - xm unique values of the dose (sorted by dose)
# - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of e
  e <- ym[1]

  # initial value of c
  # since asymptote always to the right, calculate c based on a and the max/min
  if (ad.dir == TRUE) {
    c <- max(ym) / e + 0.001 * (max(ym) / e)
  } else {
    c <- min(ym) / e - 0.001 * (min(ym) / e)
  }

  # initial value of b
  yreg <- log((ym - e * c) / (e - e * c))
  reg <- lm(yreg ~ xm)

  b <- abs(coef(reg)[2])

  startval <- list(e = e, b = b, c = c)
}


#### Exp5 ####
formExp5 <- as.formula(signal ~ e * (c - (c - 1) * exp((-1) * (b * dose)^d)))

# get starting values
startvalExp5 <- function(xm, ym, ad.dir)
# inputs
# - xm unique values of the dose (sorted by dose)
# - ym means of the signal at each value of xm (sorted by dose)
{
  # initial value of e
  e <- ym[1]

  # initial value of c
  if(ad.dir){
    c <- max(ym) / e + 0.001 * (max(ym) / e)
  }else{
    c <- min(ym) / e - 0.001 * (min(ym) / e)
  }

  # initial value of b and d
  yreg <- log((ym - e * c) / (e - e * c))
  reg <- lm(yreg ~ xm)

  b <- abs(coef(reg)[2])
  d <- exp(coef(reg)[1])

  startval <- list(e = e, b = b, c = c, d = d)
}


#### power ####
formPow <- as.formula(signal ~ e + b * (dose^c))

# get starting values
# Power model has trouble converging, so we run nls() twice
startvalPow <- function(xm, ym, ad.dir, dset) {

  require(dplyr)

  if (ad.dir) {

    e <- min(ym) - 0.001 * min(ym)

    yreg <- log(ym[xm != 0] - e)[-1]
    xreg <- log(xm[xm != 0])[-1]

    reg <- lm(yreg ~ xreg)

    c <- max(c(coef(reg)[2], 1))
    b <- exp(coef(reg)[1])

  } else {

    e <- max(ym) + 0.001 * max(ym)

    yreg <- log((ym[xm != 0] - e) * (-1))
    xreg <- log(xm[xm != 0])

    reg <- lm(yreg ~ xreg)

    c <- max(c(coef(reg)[2], 1))
    b <- exp(coef(reg)[1]) * (-1)
  }

  start_1 <- list(e = e, b = b, c = c)

  Pow <- suppressWarnings(try(nls(formPow, start = start_1, data = dset,
                                  lower = c(1, -Inf, 0.999),
                                  control = nls.control(maxiter = 500),
                                  upper = c(Inf, Inf, 18), algorithm = "port",
                                  weights = 1/(signal^2)), silent = TRUE))

  if (!inherits(Pow, "try-error")) {
    startval <- as.list(coef(Pow))
  } else {
    startval <- list(e = e, b = b, c = c)
  }

}

#### function to get bmd results
bmdres <- function(fit) {

  bmd_ci <- suppressMessages(confint(fit, "bmd"))
  bmd_mean <- coef(fit)[1] # get bmd estimate from fit
  c(bmd_mean, bmd_ci)

}


#### I use the pureErrorAnova function from alr3.
#### alr3 is now deprecated, so extracted these
#### lines from the alr3 source code in the CRAN archive. Functions were
#### minimally modified to not use S3 methods.

pureErrorAnova <- function(mod) {
  if (inherits(mod, "lm")) {
    if (is.null(mod$model)) mod <- update(mod, model = TRUE)
    p <- dim(mod$model)[2] - 1
    mod$model$Lack.of.Fit <-
      factor(randomLinComb(model.matrix(mod), 101319853))
    aov1 <- anova(mod)
    if (length(levels(mod$model$Lack.of.Fit)) == length(mod$model$Lack.of.Fit)){
      aov1
    } else {
      aov2 <- anova(
        lm(mod$model[, 1] ~ mod$model$Lack.of.Fit, weights = weights(mod))
      )
      rrow <- dim(aov1)[1]
      aov2[1, 1] <- aov1[rrow, 1] - aov2[2, 1]
      aov2[1, 2] <- aov1[rrow, 2] - aov2[2, 2]
      aov2[1, 3] <- aov2[1, 2] / aov2[1, 1]
      aov1[1:(rrow - 1), 4] <- aov1[1:(rrow - 1), 3] / aov2[2, 3]
      aov2[1, 4] <- aov2[1, 3] / aov2[2, 3]
      row.names(aov2) <- c(" Lack of fit", " Pure Error")
      aov <- rbind(aov1, aov2)
      aov[, 5] <- pf(aov[, 4], aov[, 1], aov2[2, 1], lower.tail = FALSE)
      aov
    }
  } else {
    stop("The input must be an lm object.")
  }
}

randomLinComb <- function(x, seed = NULL) {
  if (is.matrix(x)) {
    randomLinComb_matrix(x, seed)
  } else if (inherits(x, "lm")) {
    randomLinComb_lm(x, seed)
  } else {
    stop("The input must be a matrix or an lm object.")
  }
}

randomLinComb_matrix <- function(x, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  std <- function(x) {
    s <- sd(x)
    if (s > 0) (x - mean(x)) / s else x
  }
  as.vector(apply(x, 2, std) %*% as.vector(2 * rnorm(dim(x)[2]) - 1))
}

randomLinComb_lm <- function(x, seed = NULL) {
  if (is.null(x$model)) x <- update(x, model = TRUE)
  randomLinComb_matrix(x$model[, -1], seed = seed)
}


# Functions for plotting curve fits
Exp2 <- function(b, c, d, e, f, dose) {
  return(e * exp(b * dose))
}

Exp3 <- function(b, c, d, e, f, dose) {
  return(e * (exp(sign(b) * (abs(b) * dose)^d)))
}

Exp4 <- function(b, c, d, e, f, dose) {
  return(e * (c - (c - 1) * exp((-1) * b * dose)))
}

Exp5 <- function(b, c, d, e, f, dose) {
  return(e * (c - (c - 1) * exp((-1) * (b * dose)^d)))
}

Hill <- function(b, c, d, e, f, dose) {
  return(c + (d - c) / (1 + (dose / e)^b))
}

Pow <- function(b, c, d, e, f, dose) {
  return(e + b * (dose^c))
}

Poly2 <- function(b, c, d, e, f, dose) {
  return(b + c * (dose) + d * (dose^2))
}

Lin <- function(b, c, d, e, f, dose) {
  return(d + b * (dose))
}