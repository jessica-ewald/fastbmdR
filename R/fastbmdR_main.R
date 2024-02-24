# Curve fitting functions extracted from code base of FastBMD (www.fastbmd.ca)
# Please cite: https://doi.org/10.1093/bioinformatics/btaa700
# Jessica Ewald
# March 9, 2021

### NOTE: structure of "PerformCurveFitting" was based on source code of the R package DROmics.
### Specific model fits (other than Hill model); BMD calculation function; and modifications
### made throughout to be consistent with NTP approach for dose-response modeling, but
### consider citing: https://pubs.acs.org/doi/10.1021/acs.est.8b04752

### Parts of the statistical approach from here: https://doi.org/10.1371/journal.pone.0146021

### fit up to 10 different models to each gene
PerformCurveFitting <- function(data, dose,
                                ncpus = 1, models = c("Exp2","Exp3","Exp4","Exp5","Lin","Poly2","Hill","Power"))
{

  model.choices <- c("Exp2","Exp3","Exp4","Exp5","Lin","Poly2","Poly3","Poly4","Hill","Power")

  if (sum(models %in% model.choices) != length(models))
    stop("You must identify models with the correct identifiers")
  if (sum(duplicated(model.choices)) > 0)
    stop("Do not add duplicate model choices")

  requireNamespace(data.table)
  requireNamespace(dplyr)
  requireNamespace(drc)

  # definition of necessary data
  data <- as.matrix(data)
  dose <- as.numeric(dose)
  doseranks <- as.numeric(as.factor(dose))

  # get mean value for each gene, by dose
  tdata <- t(data)
  calcmean <- function(i){tapply(tdata[, i], dose, mean)}
  s <- sapply(1:nrow(data), calcmean)
  data.mean <- as.matrix(t(s))

  # calculations for starting values and other uses
  dosemin <- min(dose)
  dosemax <- max(dose)
  dosemed <- median(dose[dose!=0])
  doseu <- as.numeric(colnames(data.mean)) # sorted unique doses

  # number of points per dose-response curve
  nptsperDR <- ncol(data)
  nselect <- nrow(data)

  AICdigits <- 2 # number of digits for rounding the AIC values

  kcrit = 2 # for defining AIC

  # function to fit all the models and choose the best one
  ################################################################
  fitoneitem <- function(i)
  {
    keeplin <- "Lin" %in% models
    keepHill <- "Hill" %in% models
    keepExp2 <- "Exp2" %in% models
    keepExp3 <- "Exp3" %in% models
    keepExp4 <- "Exp4" %in% models
    keepExp5 <- "Exp5" %in% models
    keepPoly2 <- "Poly2" %in% models
    keepPoly3 <- "Poly3" %in% models
    keepPoly4 <- "Poly4" %in% models
    keepPow <- "Power" %in% models

    signal <- data[i, ]
    gene.id <- rownames(data)[i]
    signalm <- as.vector(data.mean[i,]) # means per dose
    dose0 <- signalm[1]

    # preparation of data for modelling with nls
    dset <- data.frame(dose = dose, signal = signal)
    dset <<- dset # I have to do this, otherwise the neill.test function can't find it

    # for choice of the linear trend (decreasing or increasing)
    modlin <- lm(signal ~ doseranks)
    adv.incr <- coef(modlin)[2] >= 0

    # initialize results dataframe
    item.fitres <- data.frame()

    ################## Exp2 fit ##########################
    if (keepExp2)
    {
      startExp2 <- startvalExp2(xm = doseu, ym = signalm)
      Exp2 <- suppressWarnings(try(stats::nls(formExp2, start = startExp2, data = dset,
                                       lower = c(0, -Inf), algorithm = "port"), silent = TRUE))
      if (!inherits(Exp2, "try-error"))
      {
        fit <- Exp2

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
        lof.pval.i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(Exp2)
        b.i <- par[2]
        c.i <- NA
        d.i <- NA
        e.i <- par[1]
        f.i <- NA
        SDres.i <- sigma(fit)
        mod.name <- "Exp2"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                               e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                               lof.p = lof.pval.i, ctrl.mod = ctrl)
        item.fitres <- rbind(item.fitres, res.temp)
      }
    }

    ################## Exp3 fit ##########################
    if (keepExp3)
    {
      startExp3 <- startvalExp3(xm = doseu, ym = signalm)
      Exp3 <- suppressWarnings(try(stats::nls(formExp3, start = startExp3, data = dset,
                                       lower = c(0, -Inf, -Inf), algorithm = "port"), silent = TRUE))

      if (!inherits(Exp3, "try-error"))
      {

        fit <- Exp3

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
        lof.pval.i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b.i <- par[2]
        c.i <- NA
        d.i <- par[3]
        e.i <- par[1]
        f.i <- NA
        SDres.i <- sigma(fit)
        mod.name <- "Exp3"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                               e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                               lof.p = lof.pval.i, ctrl.mod = ctrl)
        item.fitres <- rbind(item.fitres, res.temp)

      }
    }

    ################## Exp4 fit ##########################
    if (keepExp4)
    {
      startExp4 <- startvalExp4(xm = doseu, ym = signalm, ad.dir = adv.incr)

      if (adv.incr)
      {
        Exp4 <- suppressWarnings(try(stats::nls(formExp4, start = startExp4, data = dset,
                                         lower = c(1, 0, 1), algorithm = "port"), silent = TRUE))
      } else
      {
        Exp4 <- suppressWarnings(try(stats::nls(formExp4, start = startExp4, data = dset,
                                         lower = c(1, 0, 0),
                                         upper = c(Inf, Inf, 1), algorithm = "port"), silent = TRUE))
      }

      if (!inherits(Exp4, "try-error"))
      {

        fit <- Exp4

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
        lof.pval.i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b.i <- par[2]
        c.i <- par[3]
        d.i <- NA
        e.i <- par[1]
        f.i <- NA
        SDres.i <- sigma(fit)
        mod.name <- "Exp4"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                               e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                               lof.p = lof.pval.i, ctrl.mod = ctrl)
        item.fitres <- rbind(item.fitres, res.temp)

      }
    }

    ################## Exp5 fit ##########################
    if (keepExp5)
    {
      startExp5 <- startvalExp5(xm = doseu, ym = signalm, ad.dir = adv.incr)

      if (adv.incr)
      {
        Exp5 <- suppressWarnings(try(stats::nls(formExp5, start = startExp5, data = dset,
                                         lower = c(1, 0, 1, 1), algorithm = "port"), silent = TRUE))
      } else
      {
        Exp5 <- suppressWarnings(try(stats::nls(formExp5, start = startExp5, data = dset,
                                         lower = c(1, 0, 0, 1),
                                         upper = c(Inf, Inf, 1, Inf), algorithm = "port"), silent = TRUE))
      }

      if (!inherits(Exp5, "try-error"))
      {

        fit <- Exp5

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
        lof.pval.i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b.i <- par[2]
        c.i <- par[3]
        d.i <- par[4]
        e.i <- par[1]
        f.i <- NA
        SDres.i <- sigma(fit)
        mod.name <- "Exp5"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                               e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                               lof.p = lof.pval.i, ctrl.mod = ctrl)
        item.fitres <- rbind(item.fitres, res.temp)

      }
    }

    ################## Power fit ##########################
    if (keepPow)
    {
      startPow <- startvalPow(xm = doseu, ym = signalm, ad.dir = adv.incr, dset = dset)

      Pow <- suppressWarnings(try(stats::nls(formPow, start = startPow, data = dset,
                                      lower = c(1, -Inf, 0.999),
                                      upper = c(Inf, Inf, 18), algorithm = "port"), silent = TRUE))

      if (!inherits(Pow, "try-error"))
      {

        fit <- Pow

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
        lof.pval.i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b.i <- par[2]
        c.i <- par[3]
        d.i <- NA
        e.i <- par[1]
        f.i <- NA
        SDres.i <- sigma(fit)
        mod.name <- "Power"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                               e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                               lof.p = lof.pval.i, ctrl.mod = ctrl)
        item.fitres <- rbind(item.fitres, res.temp)

      }
    }

    ################## Hill fit ##########################
    if (keepHill)
    {
      startHill <- startvalHillnls2(x = dose, y = signal, xm = doseu, ym = signalm,
                                    increase = adv.incr)
      Hill <- suppressWarnings(try(stats::nls(formHill, start = startHill, data = dset,
                                       lower = c(0, -Inf, -Inf, 0), algorithm = "port"), silent = TRUE))
      if (!inherits(Hill, "try-error"))
      {

        fit <- Hill

        # collect parameters
        AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
        lof.pval.i <- neill.test(fit, dset$dose, display = FALSE)
        par <- coef(fit)
        b.i <- par["b"]
        c.i <- par["c"]
        d.i <- par["d"]
        e.i <- par["e"]
        f.i <- NA
        SDres.i <- sigma(fit)
        mod.name <- "Hill"
        ctrl <- predict(fit)[1]

        # append results to dataframe
        res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                               e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                               lof.p = lof.pval.i, ctrl.mod = ctrl)
        item.fitres <- rbind(item.fitres, res.temp)

      }
    }

    ######### Fit of the linear model ############################
    if (keeplin)
    {
      Lin <- lm(signal ~ dose, data = dset)
      fit <- Lin

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
      lof.pval.i <- pureErrorAnova(fit)[3,5]
      par <- coef(fit)
      b.i <- par[2]
      c.i <- NA
      d.i <- par[1]
      e.i <- NA
      f.i <- NA
      SDres.i <- sigma(fit)
      mod.name <- "Lin"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                             e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                             lof.p = lof.pval.i, ctrl.mod = ctrl)
      item.fitres <- rbind(item.fitres, res.temp)

    }

    ######### Fit of the Poly2 model ############################
    if (keepPoly2)
    {
      Poly2 <- lm(signal ~ dose + I(dose^2), data = dset)
      fit <- Poly2

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
      lof.pval.i <- pureErrorAnova(fit)[4,5]
      par <- coef(fit)
      b.i <- par[1]
      c.i <- par[2]
      d.i <- par[3]
      e.i <- NA
      f.i <- NA
      SDres.i <- sigma(fit)
      mod.name <- "Poly2"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                             e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                             lof.p = lof.pval.i, ctrl.mod = ctrl)
      item.fitres <- rbind(item.fitres, res.temp)

    }

    ######### Fit of the Poly3 model ############################
    if (keepPoly3)
    {
      Poly3 <- lm(signal ~ dose + I(dose^2) + I(dose^3), data = dset)
      fit <- Poly3

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
      lof.pval.i <- pureErrorAnova(fit)[5,5]
      par <- coef(fit)
      b.i <- par[1]
      c.i <- par[2]
      d.i <- par[3]
      e.i <- par[4]
      f.i <- NA
      SDres.i <- sigma(fit)
      mod.name <- "Poly3"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                             e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                             lof.p = lof.pval.i, ctrl.mod = ctrl)
      item.fitres <- rbind(item.fitres, res.temp)

    }

    ######### Fit of the Poly4 model ############################
    if (keepPoly4)
    {

      Poly4 <- lm(signal ~ dose + I(dose^2) + I(dose^3) + I(dose^4), data = dset)
      fit <- Poly4

      # collect parameters
      AIC.i <- round(AIC(fit, k = kcrit), digits = AICdigits)
      lof.pval.i <- pureErrorAnova(fit)[6,5]
      par <- coef(fit)
      b.i <- par[1]
      c.i <- par[2]
      d.i <- par[3]
      e.i <- par[4]
      f.i <- par[5]
      SDres.i <- sigma(fit)
      mod.name <- "Poly4"
      ctrl <- predict(fit)[1]

      # append results to dataframe
      res.temp <- data.frame(gene.id = gene.id, mod.name = mod.name, b = b.i, c = c.i, d = d.i,
                             e = e.i, f = f.i, SDres = SDres.i, AIC.model = AIC.i,
                             lof.p = lof.pval.i, ctrl.mod = ctrl)
      item.fitres <- rbind(item.fitres, res.temp)

    }

    item.fitres$item.ind <- c(rep(i, dim(item.fitres)[1]))
    item.fitres$ctrl.mean <- c(rep(dose0, dim(item.fitres)[1]))
    item.fitres$adv.incr <- c(rep(adv.incr, dim(item.fitres)[1]))

    rownames(item.fitres) <- NULL

    return(item.fitres)

  } ##################################### end of fitoneitem

  # Loop on items
  # parallel or sequential computation
  if (ncpus != 1) {
    clus <- parallel::makeCluster(ncpus, type = "FORK")
    res <- parallel::parLapply(clus, 1:nselect, fitoneitem)
    parallel::stopCluster(clus)
    res <- rbindlist(res)
  } else {
      res <- base::lapply(1:nselect, fitoneitem)
      res <- rbindlist(res)
  }

  dres <- as.data.frame(res)

  reslist <- list(fitres.all = dres, fitres.filt = data.frame(), data = data,
                  dose = dose, data.mean = data.mean)

  res <- structure(reslist, class = "fitObj")

  return(res);
}

#### Select the best fit model for each gene
FilterDRFit <- function(fitObj, lof.pval = 0.1)
{
  # Checks
  if (!inherits(fitObj, "fitObj"))
    stop("Use only with 'fitObj' objects, created with the function PerformCurveFitting")

  requireNamespace(data.table)
  lof.pval <- as.numeric(lof.pval)

  # get results
  fitres.all <- data.table::as.data.table(fitObj$fitres.all)
  fitres.filt <- fitres.all

  # get best fit for each feature based on selected criteria
  fitres.filt$AIC.model <- as.numeric(as.vector(fitres.filt$AIC.model))
  fitres.filt$lof.p <- as.numeric(as.vector(fitres.filt$lof.p))
  fitres.filt <- fitres.filt[as.numeric(fitres.filt$lof.p) > lof.pval]
  fitres.filt <- fitres.filt[fitres.filt[ , .I[which.min(AIC.model)], by = item.ind]$V1]


  # remove rows that had no signficant fits
  idx <- as.numeric(fitres.filt$item.ind)
  data <- fitObj$data[idx, ]
  data.mean <- fitObj$data.mean[idx, ]

  # update fit object
  fitObj$fitres.filt <- as.data.frame(fitres.filt)
  fitObj$drcfit.obj$data <- data
  fitObj$drcfit.obj$data.mean <- data.mean

  res <- structure(fitObj, class = "fitObjFilt")

  return(res)
}


### Calculation of BMD values from fitted dose-response curves
PerformBMDCalc <- function(fitObj, ncpus = 1, num.sds = 1, sample.mean = TRUE)
{

  num.sds <- as.numeric(num.sds);
  # Checks
  if (!inherits(fitObj, "fitObjFilt"))
    stop("Use only with 'fitObjFilt' objects, created with the function FilerDRFit")

  requireNamespace(data.table)
  requireNamespace(dplyr)

  dfitall <- fitObj$fitres.filt # filter this based on constant model
  dfitall$mod.name <- as.character(dfitall$mod.name)
  nselect <- nrow(dfitall)
  PerformBMDCalc
  # get necessary data for fitting
  dose <- fitObj$dose
  doseranks <- as.numeric(as.factor(dose))
  data <- fitObj$data
  data.mean <- fitObj$data.mean

  # get only correct rows
  inx.bmd <- rownames(data) %in% as.character(dfitall$gene.id)
  if(sum(inx.bmd) == 1){
    data.temp <- matrix(data[inx.bmd,], nrow = 1)
    colnames(data.temp) <- colnames(data)
    rownames(data.temp) <- rownames(data)[inx.bmd]
    data <- data.temp

    data.mean.temp <- matrix(data.mean[inx.bmd, ], nrow = 1)
    colnames(data.mean.temp) <- colnames(data.mean)
    rownames(data.mean.temp) <- rownames(data.mean)[inx.bmd]
    data.mean <- data.mean.temp
  } else {
    data <- data[inx.bmd, ]
    data.mean <- data.mean[inx.bmd, ]
  }

  item <- fitObj$fitres.filt[,1]
  fitres.bmd <- fitObj$fitres.filt

  # define BMR using either mean of controls or model evaluated at zero
  if(sample.mean){
    bmr.mode <- "ctrl.mean"
  } else {
    bmr.mode <- "ctrl.mod"
  }

  # function to calculate the bmd from the best fit model
  ################################################################
  bmdoneitem <- function(i)
  {
    # determine if adverse direction is increasing or decreasing
    adv.incr <- dfitall[i, "adv.incr"];

    # compute the BMR
    if (adv.incr == TRUE){
      bmr <- as.vector(dfitall[i, bmr.mode]) + num.sds*as.vector(dfitall[i, "SDres"])
    } else {
      bmr <- as.vector(dfitall[i, bmr.mode]) - num.sds*as.vector(dfitall[i, "SDres"])
    }

    # get best fit model type
    mod.name <- as.character(dfitall[i, "mod.name"])

    # get data for fitting
    signal <- data[i, ]
    dset <- data.frame(signal = signal, dose = dose)

    # get parameters
    b <- as.numeric(as.vector(dfitall[i, "b"]))
    c <- as.numeric(as.vector(dfitall[i, "c"]))
    d <- as.numeric(as.vector(dfitall[i, "d"]))
    e <- as.numeric(as.vector(dfitall[i, "e"]))
    f <- as.numeric(as.vector(dfitall[i, "f"]))

    # fit re-parametrized model
    switch(mod.name,
           Lin = {
             bmd0 <- (bmr - d)/b
             fit <- suppressWarnings(try(stats::nls(signal ~ d + ((bmr - d)/bmd)*dose,
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Hill = {
             bmd0 <- e*((((d-c)/(bmr-c))-1)^(1/b))
             fit <- suppressWarnings(try(stats::nls(signal ~ c + (((bmr - c)*(1 + ((bmd/e)^b)) + c) - c) / (1 + (dose/e)^b),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Exp2 = {
             bmd0 <- log(bmr/e)/b
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr/exp(b*bmd))*exp(b*dose),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Exp3 = {
             bmd0 <- ((log(bmr/e)/sign(b))^(1/d))/abs(b)
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr/exp(sign(b)*(abs(b)*bmd)^d))*(exp(sign(b)*(abs(b)*dose)^d)),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Exp4 = {
             bmd0 <- log((c-(bmr/e))/(c-1))/(-b)
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr/(c-(c-1)*exp(-1*b*bmd)))*(c - (c-1)*exp((-1)*b*dose)),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Exp5 = {
             bmd0 <- ((-log((c-(bmr/e))/(c-1)))^(1/d))/b
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr/(c-(c-1)*exp((-1)*(b*bmd)^d)))*(c - (c-1)*exp((-1)*(b*dose)^d)),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Poly2 = {
             bmd1 <- -((c + ((c^2) - 4*(b-bmr)*d)^(1/2))/(2*d))
             bmd2 <- -((c - ((c^2) - 4*(b-bmr)*d)^(1/2))/(2*d))
             bmds <- c(bmd1, bmd2)
             bmds <- bmds[bmds > 0]
             bmd0 <- min(bmds)
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr - (c*bmd + d*(bmd^2))) + c*dose + d*(dose^2),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },

           Poly3 = {
             roots <- polyroot(c(b-bmr, c, d, e))
             roots <- Re(roots)[round(Im(roots), 2) == 0]
             bmd0 <- min(roots[roots > 0])
             if(adv.incr) {bmd0 <- bmd0 - 0.1} else {bmd0 <- bmd0 + 0.1}
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr - (c*bmd + d*(bmd^2) + e*(bmd^3))) + c*dose + d*(dose^2) + e*(dose^3),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Poly4 = {
             roots <- polyroot(c(b-bmr, c, d, e, f))
             roots <- Re(roots)[round(Im(roots), 2) == 0]
             bmd0 <- min(roots[roots > 0])
             if(adv.incr) {bmd0 <- bmd0 - 0.1} else {bmd0 <- bmd0 + 0.1}
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr - (c*bmd + d*(bmd^2) + e*(bmd^3) + f*(bmd^4))) + c*dose + d*(dose^2) + e*(dose^3) + f*(dose^4),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           Power = {
             bmd0 <- ((bmr-e)/b)^(1/c)
             fit <- suppressWarnings(try(stats::nls(signal ~ (bmr - b*(bmd^c)) + b*(dose^c),
                                             data = dset, start = list(bmd = bmd0), lower = c(0), algorithm = "port"),
                                         silent = TRUE))
           },
           stop("No match here!")
    )


    # get bmd results from re-parametrized fit
    if (!inherits(fit, "try-error"))
    {
      bmd.res <- try(bmdres(fit), silent = TRUE)
    } else
    {
      bmd.res <- c(NA, NA, NA)
    }

    # replace where bmd.res failed with error code
    if ( inherits(bmd.res, "try-error")){
      bmd.res <- c(9999, NA, NA)
    }

    bmd.res <- c(bmd.res, mod.name, bmr)
    names(bmd.res) <- c("bmd", "bmdl", "bmdu", "mod.name", "bmr")
    return(bmd.res)

  } ##################################### end of bmdoneitem

  # Loop on items
  # parallel or sequential computation
  if (ncpus != 1) {
    clus <- parallel::makeCluster(ncpus, type = "FORK")
    res <- parallel::parSapply(clus, 1:nselect, bmdoneitem)
    parallel::stopCluster(clus)
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

  # change order of columns
  dres <- dres[, c("id", "mod.name", "bmd", "bmdl", "bmdu", "bmr")]

  # does bmdcalc converge for bmd, bmdl, and bmdu?
  dres$conv.pass <- rowSums(is.na(dres)) == 0

  # is the bmd < highest dose?
  dres$hd.pass <- dres$bmd < tail(dose, n = 1)

  # is the CI of the bmd narrow enough?
  dres$CI.pass <- dres$bmdu/dres$bmdl < 40

  # is the bmd < lowest dose/10?
  dres$ld.pass <- dres$bmdl > (unique(dose)[2]/10)

  # aggregate all filters
  # flag genes that don't pass low dose condition by keeping the column, but do
  # not use for the final filtering
  dres$all.pass <- (dres$conv.pass & dres$hd.pass & dres$CI.pass)

  # create output file
  dnld.file <- merge(fitObj$fitres.filt[,-12], dres[,-2], by.x = "gene.id", by.y = "id");

  return(dnld.file)

}

#### Compute POD mode ####
modePOD <- function(bmds){

  # get density plot
  density.bmd <- density(bmds, na.rm = TRUE)

  # interpolate with high # of data points
  X <- density.bmd$x
  Y <- density.bmd$y

  # calculate first derivative
  dY <- diff(Y)/diff(X)
  dX <- rowMeans(embed(X,2))

  # detect index of local maxima
  dY.signs <- sign(dY)
  dY.signs.1 <- c(0,dY.signs[-length(dY.signs)])

  # get indices of local mins and maxes
  sign.changes <- dY.signs - dY.signs.1
  inds.maxes <- which(sign.changes == -2)
  inds.mins <- which(sign.changes == 2)

  bmd.temp <- bmds
  bmds.tot <- length(bmd.temp)

  # get mins and maxes
  if((length(inds.mins) > 0) & (length(inds.maxes) > 0)){

    mins <- dX[inds.mins]
    names(mins) <- paste0("min", c(1:length(mins)))
    maxes <- dX[inds.maxes]
    names(maxes) <- paste0("max", c(1:length(maxes)))

    # order mins and maxes
    mins.maxes <- sort(c(mins,maxes))
    num.modes <- length(mins.maxes) %/% 2
    mins.maxes <- mins.maxes[1:(2*num.modes)]

    # initialize inputs to for loop
    size.peaks <- c(rep(NA, num.modes))
    per.peaks <- c(rep(NA, num.modes))

    # calculate the size of each peak
    for(i in c(1:num.modes)){

      size.peaks[i] <- sum(bmd.temp < mins.maxes[2*i])
      per.peaks[i] <- size.peaks[i]/bmds.tot

      # update bmd.temp to remove features from previous peak
      bmd.temp <- bmd.temp[!(bmd.temp < mins.maxes[2*i])]
    }

    # check size
    per.pass <- which(per.peaks > 0.05)

  } else if((length(inds.mins) == 0) & (length(inds.maxes) > 0)){

    maxes <- dX[inds.maxes]
    per.pass <- 1

  } else {
    per.pass <- NULL
  }

  # find transcriptomic pod
  if (length(per.pass > 0)){
    mode.pod <- unname(maxes[per.pass[1]])
  } else {
    mode.pod <- NA
  }
  return(mode.pod)
}

#### Compute LCRD POD ####
# as described in Crizer et al, not with Reardon modification
lcrdPOD <- function(bmds){ ## need to turn this into the actual function
  bmd.temp <- data.frame(bmd = bmds)
  bmd.temp$rank <- rank(bmd.temp$bmd)

  bmd.temp <- bmd.temp[order(bmd.temp$rank),]
  bmd.temp$ratio <- NA
  bmd.temp$ratio[1:(length(bmd.temp$ratio)-1)] <- bmd.temp$bmd[2:length(bmd.temp$ratio)]/bmd.temp$bmd[1:(length(bmd.temp$ratio)-1)]
  outlier.bmd <- which(bmd.temp$ratio > 1.66)
  if(length(outlier.bmd) == 0){
    lcrd.pod <- bmd.temp$bmd[1]
  } else {
    lcrd.pod <- bmd.temp$bmd[outlier.bmd[length(outlier.bmd)] + 1]
  }
  return(lcrd.pod)
}


#### Compute percentile-based POD ####
centilePOD <- function(bmds, centile){

  # check inputs
  stopifnot("The 'centile' parameter must be a number between 0 and 1." =
              (is.numeric(centile) & centile >= 0 & centile <= 1))

  stopifnot("The 'bmds' input must be a numeric vector with length > 1." =
              (is.numeric(bmds) & length(bmds) > 1))

  # compute pod
  cent.pod <- unname(quantile(bmds, c(centile)))

  return(cent.pod)
}


#### Compute gene set POD
gsPOD <- function(gs.lib, universe, bmd.list, overlap_number, overlap_percent) {
  # gs.lib is a list of gene sets
  # universe are all measured gene IDs
  # bmd.list is a named list of gene-level benchmark doses, with gene IDs as names
  # The gene IDs for all three must be the same type (ie. Entrez ID)

  requireNamespace(fgsea)

  # hypergeometic test
  results <- fora(gs.lib, names(bmd.list), universe) %>% suppressWarnings()
  results <- results[(size >= overlap_number) & (overlap/size >= overlap_percent)]

  if(dim(results)[1] > 0) {
    # get BMDs from each enriched gene set
    pw.bmds <- lapply(results$overlapGenes, function(x){bmd.list[unlist(x)]})
    names(pw.bmds) <- results$pathway

    # median BMD for each gene set
    bmd.med <- lapply(pw.bmds, median)
    results$bmd.med <- unlist(bmd.med)

    # pathway POD is lowest median BMD that passes all criteria
    return(min(results$bmd.med))
  } else {
    return(NA)
  }
}


rsPOD <- function(universe, bmd.list, n.iter){
  # set size = 1% size of universe
  # number sets = 1000
  # overlap must be > 5%

  res = replicate(n.iter, gsPOD(
    gs.lib = replicate(1000, sample(universe, ceiling(0.01*length(universe)), replace = F)) %>% as.data.frame() %>% as.list(),
    universe = universe,
    bmd.list = bmd.list,
    overlap_number = 0,
    overlap_percent = 0.05
  ))
  res.ci <- quantile(res, c(0.025, 0.975), na.rm = TRUE)
  res.ci <- c(POD = median(res, na.rm = T), res.ci, na_percent = signif(sum(is.na(res)/n.iter), 2))

  return(res.ci)

}

#### Compute global Mahalanobis distance ####
# As described by Nyffeler et al. (GET SPECIFIC CITATION)

prepMahalanobisDistances <- function(dat, coverVariance, treatment){
  # dat should be a matrix with samples in rows and features in columns
  # coverVariance is the percentage of variance that the retained PCs should cover
  # treatment is a vector of sample treatments, corresponding to each row in dat

  requireNamespace(tidyverse)

  ############## 1. Calculate the Eigen features from the well-level data
  PCA = prcomp(dat, center=T, scale=T) # with scaling/centering - 22 components; without scaling - 12; without both - 11
  RotationMatrix = PCA$rotation
  CumProportion = cumsum(PCA$sdev^2)/sum(PCA$sdev^2)

  ##############  2. Find the inverse of the covariance matrix

  ## Model the data to
  PC = length(which(CumProportion<coverVariance))+1
  if(PC > dim(dat)[1]){ PC = dim(dat)[1] }
  RotationMatrix <- RotationMatrix[,1:PC]
  Model = lm(PCA$x[,1:PC] ~ 0 + treatment)

  ## 3.b) get covariance matrix
  Cov = estVar(Model)

  ## 3.c) inverse
  invCov = solve(Cov)

  return(list(RotationMatrix = RotationMatrix, invCov = invCov))
}

computeMahalanobisDistance <- function(dat, RotationMatrix, invCov, treatment, controlID){
  # dat should be a matrix with samples in rows and features in columns
  # RotationMatrix are the PC loadings for the number of PCs that explain coverVariance amount of variability
  # invCov is the inverse covariance matrix computed from the PCs ~ treatment model
  # treatment is a vector with treatment labels for each row in dat
  # controlID is the treatment group that should be considered controls

  # get PC scores using loadings and number PCs computed previously
  dat = dat %*% RotationMatrix

  # compute the centroid of control samples
  ctrMean = apply(dat[treatment == controlID, ], 2, mean)

  # subtract the control centroid from each sample
  Delta = sweep(dat, 2, as.matrix(ctrMean), "-")

  # compute the Mahalanobis distance
  D = apply(Delta, 1, function(x) (x %*% invCov %*% x) %>% sqrt() %>% round(., 3))

  return(D)

}

#### Compute PLS-DA scores ####
computePLS <- function(dat, dose, ncomp){
  requireNamespace(mdatools)
  requireNamespace(dplyr)

  pls.res <- suppressWarnings(mdatools::plsda(dat, as.character(dose), ncomp = ncomp, cv = NULL, scale = T))

  # extract loadings
  xload <- pls.res$xloadings %>% as.data.frame()

  # extract scores and shift above zero
  comps <- pls.res$res$cal$xdecomp
  xscores <- comps$scores
  xscores <- t(xscores)
  xscores <- xscores + abs(min(xscores))*1.05 # shift above zero

  # extract variance explained
  varexp <- comps$expvar

  return(list(scores = xscores, loadings = xload, varExp = varexp))
}

#### Compute POD from dim reduction scores
scoresPOD <- function(dat, dose) {

  models = c("Exp2","Exp3","Exp4","Exp5","Poly2","Lin","Power","Hill")

  curve.res <- PerformCurveFitting(data = dat, dose = dose, ncpus = 1, models = models)
  curve.res <- FilterDRFit(curve.res, lof.pval = 0)

  bmds <- PerformBMDCalc(curve.res, ncpus = 1, num.sds = 1, sample.mean = TRUE)

  return(bmds)
}




