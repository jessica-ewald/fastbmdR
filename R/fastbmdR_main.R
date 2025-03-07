# FastBMD curve fitting functions
# Please cite: https://doi.org/10.1093/bioinformatics/btaa700
# Jessica Ewald

#' Structure of "PerformCurveFitting" based on source code of DROmics.
#' Specific model fits (other than Hill model);
#' BMD calculation function; and modifications
#' made throughout to be consistent with NTP approach for
#' dose-response modeling, but
#' consider citing: https://pubs.acs.org/doi/10.1021/acs.est.8b04752

#' Some methods from here: https://doi.org/10.1371/journal.pone.0146021

#' PerformCurveFitting
#'
#' This function fits the selected models for every feature in the input 'data'
#' dataframe and computes statistics describing the curve fit. It can be
#' executed in parallel by specifying 'ncpus'. It returns 'fitObj', a list that
#' includes the input data and the curve fit results.
#'
#' @importFrom parallel makeCluster parSapply stopCluster
#' @importFrom drc neill.test
#' @param data A `data.frame` with omics features in rows and sample in columns.
#' @param dose A `c()` containing the dose associated with each sample in data.
#' @param ncpus An `int` specifying the number of cpus to run in parallel.
#' @param models A `c()` of model names to use for curve fitting.
#' @return A `fitObj` object containing input data and curve fit results.
#' @export
PerformCurveFitting <- function(data,
                                dose,
                                ncpus = 1,
                                models = c("Exp2","Exp3","Exp4","Exp5","Lin","Poly2","Hill","Power")) {

  model_choices <- c("Exp2", "Exp3", "Exp4", "Exp5", "Lin",
                     "Poly2", "Poly3", "Poly4", "Hill", "Power")

  if (sum(models %in% model_choices) != length(models))
    stop("You must identify models with the correct identifiers")
  if (sum(duplicated(model_choices)) > 0)
    stop("Do not add duplicate model choices")

  # definition of necessary data
  data <- as.matrix(data)
  dose <- as.numeric(dose)
  doseranks <- as.numeric(as.factor(dose))

  # get mean value for each gene, by dose
  tdata <- t(data)
  calcmean <- function(i) {
    tapply(tdata[, i], dose, mean, na.rm = TRUE)
  }
  s <- sapply(1:nrow(data), calcmean)
  data_mean <- as.matrix(t(s))

  # calculations for starting values and other uses
  doseu <- as.numeric(colnames(data_mean)) # sorted unique doses

  # number of points per dose-response curve
  nselect <- nrow(data)

  aic_digits <- 2 # number of digits for rounding the AIC values

  kcrit <- 2 # for defining AIC

  # function to fit all the models and choose the best one
  ################################################################
  fitoneitem <- function(i) {
    keep_lin <- "Lin" %in% models
    keep_hill <- "Hill" %in% models
    keep_exp2 <- "Exp2" %in% models
    keep_exp3 <- "Exp3" %in% models
    keep_exp4 <- "Exp4" %in% models
    keep_exp5 <- "Exp5" %in% models
    keep_poly2 <- "Poly2" %in% models
    keep_poly3 <- "Poly3" %in% models
    keep_poly4 <- "Poly4" %in% models
    keep_pow <- "Power" %in% models

    signal <- data[i, ]
    gene.id <- rownames(data)[i]
    signalm <- as.vector(data_mean[i, ]) # means per dose
    dose0 <- signalm[1]

    # preparation of data for modelling with nls
    dset <- data.frame(dose = dose, signal = signal)
    dset <- na.omit(dset)
    dset <<- dset

    # for choice of the linear trend (decreasing or increasing)
    modlin <- lm(signal ~ doseranks)
    adv_incr <- coef(modlin)[2] >= 0

    # initialize results dataframe
    item_fitres <- data.frame()

    ################## Exp2 fit ##########################
    if (keep_exp2) {
      startExp2 <- startvalExp2(xm = doseu, ym = signalm)
      Exp2 <- suppressWarnings(try(nls(formExp2,
                                       start = startExp2,
                                       data = dset,
                                       lower = c(0, -Inf),
                                       algorithm = "port"),
                                   silent = TRUE))
      if (!inherits(Exp2, "try-error")) {
        fit <- Exp2

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
        lof_pval_i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(Exp2)
        b_i <- par[2]
        c_i <- NA
        d_i <- NA
        e_i <- par[1]
        f_i <- NA
        sd_res_i <- sigma(fit)
        mod_name <- "Exp2"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                               b = b_i, c = c_i, d = d_i,
                               e = e_i, f = f_i, SDres = sd_res_i,
                               AIC.model = AIC.i,
                               lof.p = lof_pval_i, ctrl.mod = ctrl)
        item_fitres <- rbind(item_fitres, res_temp)
      }
    }

    ################## Exp3 fit ##########################
    if (keep_exp3) {
      startExp3 <- startvalExp3(xm = doseu, ym = signalm)
      Exp3 <- suppressWarnings(try(nls(formExp3,
                                       start = startExp3,
                                       data = dset,
                                       lower = c(0, -Inf, -Inf),
                                       algorithm = "port"),
                                   silent = TRUE))

      if (!inherits(Exp3, "try-error")) {
        fit <- Exp3

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
        lof_pval_i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b_i <- par[2]
        c_i <- NA
        d_i <- par[3]
        e_i <- par[1]
        f_i <- NA
        sd_res_i <- sigma(fit)
        mod_name <- "Exp3"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                               b = b_i, c = c_i, d = d_i,
                               e = e_i, f = f_i, SDres = sd_res_i,
                               AIC.model = AIC.i,
                               lof.p = lof_pval_i, ctrl.mod = ctrl)
        item_fitres <- rbind(item_fitres, res_temp)

      }
    }
    ################## Exp4 fit ##########################
    if (keep_exp4) {
      startExp4 <- startvalExp4(xm = doseu, ym = signalm, ad.dir = adv_incr)

      if (adv_incr) {
        Exp4 <- suppressWarnings(try(nls(formExp4,
                                         start = startExp4,
                                         data = dset,
                                         lower = c(1, 0, 1),
                                         algorithm = "port"),
                                     silent = TRUE))
      } else {
        Exp4 <- suppressWarnings(try(nls(formExp4,
                                         start = startExp4,
                                         data = dset,
                                         lower = c(1, 0, 0),
                                         upper = c(Inf, Inf, 1),
                                         algorithm = "port"),
                                     silent = TRUE))
      }

      if (!inherits(Exp4, "try-error")) {
        fit <- Exp4

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
        lof_pval_i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b_i <- par[2]
        c_i <- par[3]
        d_i <- NA
        e_i <- par[1]
        f_i <- NA
        sd_res_i <- sigma(fit)
        mod_name <- "Exp4"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                               b = b_i, c = c_i, d = d_i,
                               e = e_i, f = f_i, SDres = sd_res_i,
                               AIC.model = AIC.i,
                               lof.p = lof_pval_i, ctrl.mod = ctrl)
        item_fitres <- rbind(item_fitres, res_temp)

      }
    }
    ################## Exp5 fit ##########################
    if (keep_exp5) {
      startExp5 <- startvalExp5(xm = doseu, ym = signalm, ad.dir = adv_incr)

      if (adv_incr) {
        Exp5 <- suppressWarnings(try(nls(formExp5,
                                         start = startExp5,
                                         data = dset,
                                         lower = c(1, 0, 1, 1),
                                         algorithm = "port"),
                                     silent = TRUE))
      } else {
        Exp5 <- suppressWarnings(try(nls(formExp5,
                                         start = startExp5,
                                         data = dset,
                                         lower = c(1, 0, 0, 1),
                                         upper = c(Inf, Inf, 1, Inf),
                                         algorithm = "port"),
                                     silent = TRUE))
      }

      if (!inherits(Exp5, "try-error")) {
        fit <- Exp5

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
        lof_pval_i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b_i <- par[2]
        c_i <- par[3]
        d_i <- par[4]
        e_i <- par[1]
        f_i <- NA
        sd_res_i <- sigma(fit)
        mod_name <- "Exp5"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                               b = b_i, c = c_i, d = d_i,
                               e = e_i, f = f_i, SDres = sd_res_i,
                               AIC.model = AIC.i,
                               lof.p = lof_pval_i, ctrl.mod = ctrl)
        item_fitres <- rbind(item_fitres, res_temp)
      }
    }

    ################## Power fit ##########################
    if (keep_pow) {
      startPow <- startvalPow(xm = doseu,
                              ym = signalm,
                              ad.dir = adv_incr,
                              dset = dset)

      Pow <- suppressWarnings(try(nls(formPow,
                                      start = startPow,
                                      data = dset,
                                      lower = c(1, -Inf, 0.999),
                                      upper = c(Inf, Inf, 18),
                                      algorithm = "port"),
                                  silent = TRUE))

      if (!inherits(Pow, "try-error")) {
        fit <- Pow

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
        lof_pval_i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b_i <- par[2]
        c_i <- par[3]
        d_i <- NA
        e_i <- par[1]
        f_i <- NA
        sd_res_i <- sigma(fit)
        mod_name <- "Power"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                               b = b_i, c = c_i, d = d_i,
                               e = e_i, f = f_i, SDres = sd_res_i,
                               AIC.model = AIC.i,
                               lof.p = lof_pval_i, ctrl.mod = ctrl)
        item_fitres <- rbind(item_fitres, res_temp)

      }
    }

    ################## Hill fit ##########################
    if (keep_hill) {
      startHill <- startvalHillnls2(x = dose,
                                    y = signal,
                                    xm = doseu,
                                    ym = signalm,
                                    increase = adv_incr)
      Hill <- suppressWarnings(try(nls(formHill,
                                       start = startHill,
                                       data = dset, 
                                       lower = c(0, -Inf, -Inf, 0),
                                       algorithm = "port"),
                                   silent = TRUE))
      if (!inherits(Hill, "try-error")) {
        fit <- Hill

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
        lof_pval_i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b_i <- par["b"]
        c_i <- par["c"]
        d_i <- par["d"]
        e_i <- par["e"]
        f_i <- NA
        sd_res_i <- sigma(fit)
        mod_name <- "Hill"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                               b = b_i, c = c_i, d = d_i,
                               e = e_i, f = f_i, SDres = sd_res_i,
                               AIC.model = AIC.i,
                               lof.p = lof_pval_i, ctrl.mod = ctrl)
        item_fitres <- rbind(item_fitres, res_temp)
      }
    }
    ######### Fit of the linear model ############################
    if (keep_lin) {
      Lin <- lm(signal ~ dose, data = dset)
      fit <- Lin

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
      lof_pval_i <- pureErrorAnova(fit)[3, 5]
      par <- coef(fit)
      b_i <- par[2]
      c_i <- NA
      d_i <- par[1]
      e_i <- NA
      f_i <- NA
      sd_res_i <- sigma(fit)
      mod_name <- "Lin"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      # append results to dataframe
      res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                             b = b_i, c = c_i, d = d_i,
                             e = e_i, f = f_i, SDres = sd_res_i,
                             AIC.model = AIC.i,
                             lof.p = lof_pval_i, ctrl.mod = ctrl)
      item_fitres <- rbind(item_fitres, res_temp)
    }

    ######### Fit of the Poly2 model ############################
    if (keep_poly2) {
      Poly2 <- lm(signal ~ dose + I(dose^2), data = dset)
      fit <- Poly2

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
      lof_pval_i <- pureErrorAnova(fit)[4, 5]
      par <- coef(fit)
      b_i <- par[1]
      c_i <- par[2]
      d_i <- par[3]
      e_i <- NA
      f_i <- NA
      sd_res_i <- sigma(fit)
      mod_name <- "Poly2"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      # append results to dataframe
      res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                             b = b_i, c = c_i, d = d_i,
                             e = e_i, f = f_i, SDres = sd_res_i,
                             AIC.model = AIC.i,
                             lof.p = lof_pval_i, ctrl.mod = ctrl)
      item_fitres <- rbind(item_fitres, res_temp)
    }
    ######### Fit of the Poly3 model ############################
    if (keep_poly3) {
      Poly3 <- lm(signal ~ dose + I(dose^2) + I(dose^3), data = dset)
      fit <- Poly3

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
      lof_pval_i <- pureErrorAnova(fit)[5,5]
      par <- coef(fit)
      b_i <- par[1]
      c_i <- par[2]
      d_i <- par[3]
      e_i <- par[4]
      f_i <- NA
      sd_res_i <- sigma(fit)
      mod_name <- "Poly3"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                             b = b_i, c = c_i, d = d_i,
                             e = e_i, f = f_i, SDres = sd_res_i,
                             AIC.model = AIC.i,
                             lof.p = lof_pval_i, ctrl.mod = ctrl)
      item_fitres <- rbind(item_fitres, res_temp)
    }

    ######### Fit of the Poly4 model ############################
    if (keep_poly4) {
      Poly4 <- lm(signal ~ dose + I(dose^2) + I(dose^3) + I(dose^4),
                  data = dset)
      fit <- Poly4

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = aic_digits)
      lof_pval_i <- pureErrorAnova(fit)[6, 5]
      par <- coef(fit)
      b_i <- par[1]
      c_i <- par[2]
      d_i <- par[3]
      e_i <- par[4]
      f_i <- par[5]
      sd_res_i <- sigma(fit)
      mod_name <- "Poly4"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      res_temp <- data.frame(gene.id = gene.id, mod.name = mod_name,
                             b = b_i, c = c_i, d = d_i,
                             e = e_i, f = f_i, SDres = sd_res_i,
                             AIC.model = AIC.i,
                             lof.p = lof_pval_i, ctrl.mod = ctrl)
      item_fitres <- rbind(item_fitres, res_temp)

    }
    item_fitres$item.ind <- c(rep(i, dim(item_fitres)[1]))
    item_fitres$ctrl.mean <- c(rep(dose0, dim(item_fitres)[1]))
    item_fitres$adv.incr <- c(rep(adv_incr, dim(item_fitres)[1]))

    rownames(item_fitres) <- NULL
    return(item_fitres)
  } ##################################### end of fitoneitem

  # Loop on items
  # parallel or sequential computation
  if (ncpus != 1) {
    clus <- makeCluster(ncpus, type = "FORK")
    res <- parLapply(clus, 1:nselect, fitoneitem)
    stopCluster(clus)
    res <- rbindlist(res)
  } else {
    res <- base::lapply(1:nselect, fitoneitem)
    res <- rbindlist(res)
  }

  dres <- as.data.frame(res)
  reslist <- list(fitres.all = dres, fitres.filt = data.frame(), data = data,
                  dose = dose, data_mean = data_mean)
  res <- structure(reslist, class = "fitObj")

  return(res)
}


#' FilterDRFit
#'
#' This function filters curves according to specified qaqc thresholds and
#' selects the best-fit curve from the remaining curves for each feature
#' according to the specified filt.var variable. It returns an updated fitObj
#' object that includes the filtered curve fits.
#'
#' @import data.table
#' @param fitObj A `fitObj` object created by PerformCurveFitting.
#' @param lof.pval A `numeric` object that specifies the lack-of-fit p-value threshold.
#' @param filt.var A `character` object that specifies the variable used for selecting the best curve fit.
#' @return An updated `fitObj` object with filtered curve fits.
#' @export
FilterDRFit <- function(fitObj, lof.pval = 0.1, filt.var = "AIC.model") {
  # Checks
  if (!inherits(fitObj, "fitObj"))
    stop("Use only with 'fitObj' objects, 
          created with the function PerformCurveFitting")

  lof.pval <- as.numeric(lof.pval)

  # get results
  fitres.all <- as.data.table(fitObj$fitres.all)
  fitres.filt <- fitres.all

  # get best fit for each feature based on selected criteria
  fitres.filt$AIC.model <- as.numeric(as.vector(fitres.filt$AIC.model))
  fitres.filt$SDres <- as.numeric(as.vector(fitres.filt$SDres))
  fitres.filt$lof.p <- as.numeric(as.vector(fitres.filt$lof.p))
  fitres.filt <- fitres.filt[as.numeric(fitres.filt$lof.p) > lof.pval]
  fitres.filt <- fitres.filt[fitres.filt[,
                                         .I[which.min(get(filt.var))],
                                         by = item.ind]$V1]


  # remove rows that had no signficant fits
  idx <- as.numeric(fitres.filt$item.ind)
  data <- fitObj$data[idx, ]
  data_mean <- fitObj$data_mean[idx, ]

  # update fit object
  fitObj$fitres.filt <- as.data.frame(fitres.filt)
  fitObj$drcfit.obj$data <- data
  fitObj$drcfit.obj$data_mean <- data_mean

  res <- structure(fitObj, class = "fitObjFilt")

  return(res)
}

#' PerformBMDCalc
#'
#' This function computes the benchmark dose (BMD) for each curve fit. The BMD
#' is defined as where the fitted curve reaches the benchmark response (BMR).
#' The BMR is defined relative to the controls values, for example as the mean
#' of the control values, plus or minus 1 standard deviation of the controls.
#' Various BMD QA/QC metrics are computed to assess BMD uncertainty. This
#' function returns a dataframe with the BMD fit results.
#'
#' @importFrom parallel makeCluster parSapply stopCluster
#' @param fitObj A `fitObj` object from FilterDRFit.
#' @param ncpus An `int` specifying the number of cpus to use in parallel.
#' @param num.sds A `numeric` specifying the number of std to define the BMR.
#' @param bmr.method A `character` specifying the BMR definition method.
#' @param log10.dose A `logical` specifying whether the dose is on a log scale.
#' @return A `data.frame` with curve fit and BMD results.
#' @export
PerformBMDCalc <- function(fitObj, ncpus = 1, num.sds = 1,
                           bmr.method = "sample.mean", log10.dose = FALSE)
{

  num.sds <- as.numeric(num.sds)

  if (!inherits(fitObj, "fitObjFilt"))
    stop("Use only with 'fitObjFilt' objects, 
          created with the function FilterDRFit")

  dfitall <- fitObj$fitres.filt # filter this based on constant model
  dfitall$mod.name <- as.character(dfitall$mod.name)
  nselect <- nrow(dfitall)

  # get necessary data for fitting
  dose <- fitObj$dose
  doseranks <- as.numeric(as.factor(dose))
  data <- fitObj$data
  data_mean <- fitObj$data_mean

  # get only correct rows
  inx_bmd <- rownames(data) %in% as.character(dfitall$gene.id)
  if (sum(inx_bmd) == 1) {
    data_temp <- matrix(data[inx_bmd,], nrow = 1)
    colnames(data_temp) <- colnames(data)
    rownames(data_temp) <- rownames(data)[inx_bmd]
    data <- data_temp

    data_mean_temp <- matrix(data_mean[inx_bmd, ], nrow = 1)
    colnames(data_mean_temp) <- colnames(data_mean)
    rownames(data_mean_temp) <- rownames(data_mean)[inx_bmd]
    data_mean <- data_mean_temp
  } else {
    data <- data[inx_bmd, ]
    data_mean <- data_mean[inx_bmd, ]
  }

  item <- fitObj$fitres.filt[, 1]
  fitres.bmd <- fitObj$fitres.filt

  # define BMR using either mean of controls or model evaluated at zero
  if(bmr.method == "sample.mean") {
    bmr.mode <- "ctrl.mean"
  } else if(bmr.method == "model.mean") {
    bmr.mode <- "ctrl.mod"
  } else {
    one_sided_pval <- (1 - pnorm(abs(num.sds)))
    bmr.mode <- "ctrl.quantile"
  }

  # function to calculate the bmd from the best fit model
  ################################################################
  bmdoneitem <- function(i)
  {
    # determine if adverse direction is increasing or decreasing
    adv_incr <- dfitall[i, "adv.incr"]

    # get best fit model type
    mod.name <- as.character(dfitall[i, "mod.name"])

    # get data for fitting
    signal <- data[i, ]
    dset <- data.frame(signal = signal, dose = dose)
    dset <- na.omit(dset)
    ctrl_vals <- signal[dose == 0]
    SDctrl <- sd(ctrl_vals)

    # compute the BMR
    if (bmr.mode == "ctrl.quantile") {
      # specify quantile-based bmr
      if (adv_incr == TRUE) {
        bmr <- quantile(ctrl_vals, (1 - one_sided_pval), na.rm = TRUE)
      } else {
        bmr <- quantile(ctrl_vals, one_sided_pval, na.rm = TRUE)
      }
    } else {
      sdres <- as.vector(dfitall[i, "SDres"])
      if (adv_incr == TRUE) {
        bmr <- as.vector(dfitall[i, bmr.mode]) + num.sds * sdres
      } else {
        bmr <- as.vector(dfitall[i, bmr.mode]) - num.sds * sdres
      }
    }

    # get parameters
    b <- as.numeric(as.vector(dfitall[i, "b"]))
    c <- as.numeric(as.vector(dfitall[i, "c"]))
    d <- as.numeric(as.vector(dfitall[i, "d"]))
    e <- as.numeric(as.vector(dfitall[i, "e"]))
    f <- as.numeric(as.vector(dfitall[i, "f"]))

    # fit re-parametrized model
    switch(mod.name,
      Lin = {
        bmd0 <- (bmr - d) / b
        rf <- signal ~ d + ((bmr - d) / bmd) * dose
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Hill = {
        bmd0 <- e * ((((d - c) / (bmr - c)) - 1)^(1 / b))
        rf <- signal ~ c + (((bmr - c) * (1 + ((bmd / e)^b)) + c) - c) / (1 + (dose / e)^b)
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Exp2 = {
        bmd0 <- log(bmr / e) / b
        rf <- signal ~ (bmr / exp(b * bmd)) * exp(b * dose)
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Exp3 = {
        bmd0 <- ((log(bmr / e) / sign(b))^(1 / d)) / abs(b)
        rf <- signal ~ (bmr / exp(sign(b) * (abs(b) * bmd)^d)) * (exp(sign(b) * (abs(b) * dose)^d))
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Exp4 = {
        bmd0 <- log((c - (bmr / e)) / (c - 1)) / (-b)
        rf <- signal ~ (bmr/(c - (c - 1) * exp(-1 * b * bmd))) * (c - (c - 1) * exp((-1) * b * dose))
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Exp5 = {
        bmd0 <- ((-log((c - (bmr / e)) / (c - 1)))^(1 / d)) / b
        rf <- signal ~ (bmr / (c - (c - 1) * exp((-1) * (b * bmd)^d))) * (c - (c - 1) * exp((-1) * (b * dose)^d))
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Poly2 = {
        bmd1 <- -((c + ((c^2) - 4 * (b - bmr) * d)^(1 / 2)) / (2 * d))
        bmd2 <- -((c - ((c^2) - 4 * (b - bmr) * d)^(1 / 2)) / (2 * d))
        bmds <- c(bmd1, bmd2)
        bmds <- bmds[bmds > 0]
        bmd0 <- min(bmds)
        rf <- signal ~ (bmr - (c * bmd + d * (bmd^2))) + c * dose + d * (dose^2)
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Poly3 = {
        roots <- polyroot(c(b - bmr, c, d, e))
        roots <- Re(roots)[round(Im(roots), 2) == 0]
        bmd0 <- min(roots[roots > 0])
        if (adv_incr) {
          bmd0 <- bmd0 - 0.1
        } else {
          bmd0 <- bmd0 + 0.1
        }
        rf <- signal ~ (bmr - (c * bmd + d * (bmd^2) + e * (bmd^3))) + c * dose + d * (dose^2) + e * (dose^3)
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Poly4 = {
        roots <- polyroot(c(b - bmr, c, d, e, f))
        roots <- Re(roots)[round(Im(roots), 2) == 0]
        bmd0 <- min(roots[roots > 0])
        if (adv_incr) {
          bmd0 <- bmd0 - 0.1
        } else {
          bmd0 <- bmd0 + 0.1
        }
        rf <- signal ~ (bmr - (c * bmd + d * (bmd^2) + e * (bmd^3) + f * (bmd^4))) + c * dose + d * (dose^2) + e * (dose^3) + f * (dose^4)
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      Power = {
        bmd0 <- ((bmr - e)/b)^(1 / c)
        rf <- signal ~ (bmr - b * (bmd^c)) + b * (dose^c)
        fit <- suppressWarnings(try(nls(rf,
                                        data = dset,
                                        start = list(bmd = bmd0),
                                        lower = c(0),
                                        algorithm = "port"),
                                    silent = TRUE))
      },
      stop("No match here!")
    )


    # get bmd results from re-parametrized fit
    if (!inherits(fit, "try-error")) {
      bmd.res <- try(bmdres(fit), silent = TRUE)
    } else {
      bmd.res <- c(NA, NA, NA)
    }

    # replace where bmd.res failed with error code
    if (inherits(bmd.res, "try-error")) {
      bmd.res <- c(9999, NA, NA)
    }

    bmd.res <- c(bmd.res, mod.name, bmr, SDctrl)
    names(bmd.res) <- c("bmd", "bmdl", "bmdu", "mod.name", "bmr", "SDctrl")
    return(bmd.res)

  } ##################################### end of bmdoneitem

  # Loop on items
  # parallel or sequential computation
  if (ncpus != 1) {
    clus <- makeCluster(ncpus, type = "FORK")
    res <- parSapply(clus, 1:nselect, bmdoneitem)
    stopCluster(clus)
  } else {
    res <- sapply(1:nselect, bmdoneitem)
  }

  # change class of columns
  dres <- as.data.frame(t(res))
  dres$bmd <- as.numeric(as.character(dres$bmd))
  dres$bmdl <- as.numeric(as.character(dres$bmdl))
  dres$bmdu <- as.numeric(as.character(dres$bmdu))
  dres$id <- as.character(item)
  dres$bmr <- as.numeric(as.character(dres$bmr))
  dres$SDctrl <- as.numeric(as.character(dres$SDctrl))

  # change order of columns
  dres <- dres[, c("id", "mod.name", "bmd", "bmdl", "bmdu", "bmr", "SDctrl")]

  # does bmdcalc converge for bmd, bmdl, and bmdu?
  dres$conv.pass <- rowSums(is.na(dres)) == 0

  # is the bmd < highest dose?
  dres$hd.pass <- dres$bmd < max(dose, n = 1)

  # is the CI of the bmd narrow enough?
  # is the BMD much lower than the lowest dose?
  if (log10.dose == TRUE) {
    dres$CI.pass <- dres$bmdu - dres$bmdl < log10(40)
    dres$ld.pass <- dres$bmdl > (unique(sort(dose))[2] - log10(10))
  } else {
    dres$CI.pass <- dres$bmdu/dres$bmdl < 40
    dres$ld.pass <- dres$bmdl > (unique(sort(dose))[2] / 10)
  }

  # aggregate all filters
  # flag genes that don't pass low dose condition by keeping the column, but do
  # not use for the final filtering
  dres$all.pass <- (dres$conv.pass & dres$hd.pass & dres$CI.pass)

  # create output file
  dnld_file <- merge(fitObj$fitres.filt[, -12],
                     dres[, -2],
                     by.x = "gene.id",
                     by.y = "id")

  return(dnld_file)

}


#' scoresPOD
#'
#' This function around the FastBMD workflow with parameters set to work well
#' for omics datasets that have been summarized as Mahalanobis distances. In
#' this scenario, it is important to be extra conservative with filtering curves
#' out because there is only a small number of features (one for global, ~15 for
#' categorical) thus losing marginal fits can greatly impact the results. Also,
#' a symmetric BMD (+/- SD) is not appropriate because distance values cannot go
#' below zero.
#'
#' @param dat A `data.frame` with omics features in rows and sample in columns.
#' @param dose A `c()` containing the dose associated with each sample in data.
#' @param log10.dose A `logical` specifying whether the dose is on a log scale.
#' @param num.sds A `numeric` specifying the number of std to define the BMR.
#' @param filt.var A `character` that specifies how to select the best fit.
#' @return A `data.frame` with curve fit and BMD results.
#' @export
scoresPOD <- function(dat, dose, log10.dose = FALSE,
                      num.sds = 1, filt.var = "AIC.model") {

  models <- c("Exp2", "Exp3", "Exp4", "Exp5", "Poly2", "Lin", "Power", "Hill")

  curve_res <- PerformCurveFitting(data = dat, dose = dose,
                                   ncpus = 1, models = models)
  curve_res <- FilterDRFit(curve_res, lof.pval = 0, filt.var = filt.var)

  if (dim(curve_res$fitres.filt)[1] > 0) {
    bmds <- PerformBMDCalc(curve_res, ncpus = 1, num.sds = num.sds,
                           bmr.method = "ctrl.quantile",
                           log10.dose = log10.dose)
  } else {
    bmds <- NULL
  }

  return(bmds)
}