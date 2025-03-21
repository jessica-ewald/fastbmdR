#' plot_bmd_curve
#'
#' This function plots both the observed data points and the fitted curve for
#' one feature in the fitres.filt dataframe in a fitObj.
#'
#' @importFrom ggplot ggplot geom_point geom_line geom_vline xlab ylab xlim
#' @importFrom ggplot theme theme_bw aes element_text
#' @param feature A `character` specifying which feature to plot a curve for.
#' @param fitObj A `fitObj` containing the curve fits and data.
#' @param dose_spacing A `character` specifying the scale of the dose range.
#' @param return_type A `character` specifying the return data format.
#' @return Either a plot, ggplot grob, or dataframe.
#' @export
plot_bmd_curve <- function(feature,
                           fitObj,
                           dose_spacing = "natural",
                           return_type = "plot") {

  # Checks
  if (!inherits(fitObj, "fitObjFiltBMD"))
    stop("Use only with 'fitObjFiltBMD' objects, 
          created with the function PerformBMDCalc")

  # get relevant objects from fitObj
  bmd_df <- fitObj$bmd_res
  dat <- as.data.frame(fitObj$data)
  dose <- fitObj$dose

  # extract parameters
  plot_fit <- bmd_df[bmd_df$gene.id == feature, ]
  model <- plot_fit$mod.name
  b <- plot_fit$b
  c <- plot_fit$c
  d <- plot_fit$d
  e <- plot_fit$e
  f <- plot_fit$f
  bmd <- plot_fit$bmd
  bmdl <- plot_fit$bmdl
  bmdu <- plot_fit$bmdu

  # interpolate curve fit
  highest_dose <- max(unique(dose))
  if (dose_spacing == "natural") {

    x <- seq(0, highest_dose, length.out = 100)

  } else if (dose_spacing == "log10") {

    x <- (10^(seq(0, log10(400), length.out = 100))) - 1

  } else if (dose_spacing == "log2") {

    x <- (2^(seq(0, log2(400), length.out = 100))) - 1

  } else {

    return("dose_spacing parameter value must be either 
            'natural', 'log10', or 'log2'")

  }

  f_x <- switch(model,
                "Exp2" = Exp2(b, c, d, e, f, x),
                "Exp3" = Exp3(b, c, d, e, f, x),
                "Exp4" = Exp4(b, c, d, e, f, x),
                "Exp5" = Exp5(b, c, d, e, f, x),
                "Hill" = Hill(b, c, d, e, f, x),
                "Pow" = Pow(b, c, d, e, f, x),
                "Poly2" = Poly2(b, c, d, e, f, x),
                "Lin" = Lin(b, c, d, e, f, x))
  plot.interp <- data.frame(x = x, f_x = f_x)


  # Get observed data
  plot.dat <- as.data.frame(t(dat[feature, ]))
  plot.dat$x <- dose
  plot.dat <- plot.dat[, c("x", feature)]
  colnames(plot.dat) <- c("x", "Observations")

  # Merge interpolated and observed data
  plot.results <- merge(plot.dat, plot.interp,
                        by = "x", all = TRUE)
  # Make plot
  p <- ggplot(plot.results, aes(x = x)) +
    geom_point(aes(y = Observations)) +
    geom_line(aes(y = f_x)) +
    geom_vline(aes(xintercept = bmdl),
               linetype = "dashed",
               color = "red",
               na.rm = TRUE) +
    geom_vline(aes(xintercept = bmd),
               linetype = "solid",
               color = "red",
               na.rm = TRUE) +
    geom_vline(aes(xintercept = bmdu),
               linetype = "dashed",
               color = "red",
               na.rm = TRUE) +
    xlim(0, highest_dose) +
    xlab("Concentration") +
    ylab("Response") +
    theme_bw() +
    theme(strip.text = element_text(size = 6))

  if (return_type == "plot") {

    print(p)

  } else if (return_type == "plot.object") {

    return(p)

  } else if (return_type == "plot.data") {

    return(plot.results)

  } else {

    return("return_type parameter value must be either 
            'plot', 'plot.object', or 'plot.data'")

  }
}


# Helper functions for curve fits

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
