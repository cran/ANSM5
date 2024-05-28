#' Calculate Theil-Kendall beta
#'
#' @description
#' `theil.kendall()` calculates the Theil-Kendall beta and is used in chapter 11 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Numeric vector of same length as x
#' @param x Numeric vector of same length as y
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param do.abbreviated Boolean indicating whether or not to use abbreviated Theil procedure (defaults to `FALSE`)
#' @param do.alpha Boolean indicating whether or not to report estimate of alpha (defaults to `FALSE`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `FALSE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMstat object with the results from applying the function
#' @examples
#' # Example 11.6 from "Applied Nonparametric Statistical Methods" (5th edition)
#' theil.kendall(ch11$reportedtime, ch11$parentlimit, do.alpha = TRUE)
#'
#' # Exercise 11.10 from "Applied Nonparametric Statistical Methods" (5th edition)
#' theil.kendall(ch11$N.Scotland, ch11$SW.England)
#'
#' @importFrom stats complete.cases median pnorm quantile
#' @export
theil.kendall <-
  function(y, x, H0 = NULL, do.abbreviated = FALSE, do.alpha = FALSE,
           alternative = c("two.sided", "less", "greater"), CI.width = 0.95,
           max.exact.cases = 10, nsims.mc = 100000, seed = NULL,
           do.asymp = FALSE, do.exact = TRUE, do.CI = FALSE, do.mc = FALSE) {
    stopifnot(is.numeric(y), is.numeric(x), length(y) == length(x),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.logical(do.abbreviated) == TRUE, is.logical(do.alpha) == TRUE,
              is.numeric(CI.width), length(CI.width) == 1,
              CI.width > 0, CI.width < 1,
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE, is.logical(do.mc) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- paste0(deparse(substitute(y)), " ~ ", deparse(substitute(x)))
    varname2 <- NULL
    varname3 <- NULL

    #unused arguments
    cont.corr <- NULL
    #default outputs
    pval <- NULL
    pval.stat <- NULL
    pval.note <- NULL
    pval.asymp <- NULL
    pval.asymp.stat <- NULL
    pval.asymp.note <- NULL
    pval.exact <- NULL
    pval.exact.stat <- NULL
    pval.exact.note <- NULL
    pval.mc <- NULL
    pval.mc.stat <- NULL
    pval.mc.note <- NULL
    actualCIwidth.exact <- NULL
    CI.exact.lower <- NULL
    CI.exact.upper <- NULL
    CI.exact.note <- NULL
    CI.asymp.lower <- NULL
    CI.asymp.upper <- NULL
    CI.asymp.note <- NULL
    CI.mc.lower <- NULL
    CI.mc.upper <- NULL
    CI.mc.note <- NULL
    CI.sample.lower <- NULL
    CI.sample.upper <- NULL
    CI.sample.note <- NULL
    stat.note <- NULL

    #prepare
    complete.cases.id <- complete.cases(x, y)
    y <- y[complete.cases.id] #remove missing cases
    x <- x[complete.cases.id] #remove missing cases
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    x.c <- x - median(x)
    n <- length(y)

    #calculate estimate of beta
    outer.y <- outer(y, y, "-")
    outer.x <- outer(x.c, x.c, "-")
    bvals <- outer.y / outer.x
    if (do.abbreviated){
      bvals2 <- NULL
      if (n %% 2 == 0){ #even
        for (i in 1:(n / 2)){
          bvals2 <- c(bvals2, bvals[i, i + (n / 2)])
        }
      }else{ #odd
        for (i in 1:((n - 1) / 2)){
          bvals2 <- c(bvals2, bvals[i, i + ((n +1) / 2)])
        }
      }
    }else{
      bvals2 <- bvals[upper.tri(bvals)]
    }
    #is.finite used to ignore tied x values giving infinite slopes
    stat <- median(bvals2[is.finite(bvals2)])
    statlabel <- "Theil-Kendall beta"
    if (!is.null(H0)){
      Tt <- sum(sign(bvals2 - H0) == 1) - sum(sign(bvals2 - H0) == -1)
    }

    #give mc output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.mc <- TRUE
    }

    #exact p-value
    if(do.exact && n <= max.exact.cases && !is.null(H0)){
      outer.y <- outer(y, y, "-")
      permutations <- perms(n)
      n.perms <- dim(permutations)[1]
      pval.exact <- 0
      for (i in 1:n.perms){
        outer.x <- outer(x.c[permutations[i,]], x.c[permutations[i,]], "-")
        bvals <- outer.y / outer.x
        if (do.abbreviated){
          bvals2 <- NULL
          if (n %% 2 == 0){ #even
            for (i in 1:(n / 2)){
              bvals2 <- c(bvals2, bvals[i, i + (n / 2)])
            }
          }else{ #odd
            for (i in 1:((n - 1) / 2)){
              bvals2 <- c(bvals2, bvals[i, i + ((n +1) / 2)])
            }
          }
        }else{
          bvals2 <- bvals[upper.tri(bvals)]
        }
        #is.finite used to ignore tied x values giving infinite slopes
        beta.tmp <- median(bvals2[is.finite(bvals2)])
        if (alternative == "two.sided"){
          if (abs(beta.tmp) >= abs(stat)){
            pval.exact <- pval.exact + 1 / n.perms
          }
        }else if (alternative == "less"){
          if (beta.tmp >= stat){
            pval.exact <- pval.exact + 1 / n.perms
          }
        }else if (alternative == "greater"){
          if (beta.tmp <= stat){
            pval.exact <- pval.exact + 1 / n.perms
          }
        }
      }
    }

    #Monte Carlo p-value
    if(do.mc && !is.null(H0)){
      if (!is.null(seed)){set.seed(seed)}
      outer.y <- outer(y, y, "-")
      pval.mc <- 0
      for (i in 1:nsims.mc){
        x.tmp <- x.c[sample(n, n, replace = FALSE)]
        outer.x <- outer(x.tmp, x.tmp, "-")
        bvals <- outer.y / outer.x
        if (do.abbreviated){
          bvals2 <- NULL
          if (n %% 2 == 0){ #even
            for (i in 1:(n / 2)){
              bvals2 <- c(bvals2, bvals[i, i + (n / 2)])
            }
          }else{ #odd
            for (i in 1:((n - 1) / 2)){
              bvals2 <- c(bvals2, bvals[i, i + ((n +1) / 2)])
            }
          }
        }else{
          bvals2 <- bvals[upper.tri(bvals)]
        }
        #is.finite used to ignore tied x values giving infinite slopes
        beta.tmp <- median(bvals2[is.finite(bvals2)])
        if (alternative == "two.sided"){
          if (abs(beta.tmp) >= abs(stat)){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }else if (alternative == "less"){
          if (beta.tmp >= stat){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }else if (alternative == "greater"){
          if (beta.tmp <= stat){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }
      }
    }

    #asymptotic p-value
    if(do.asymp && !do.abbreviated && !is.null(H0)){
      var <- n * (n - 1) * (2 * n + 5) / 18
      z <- Tt / sqrt(var)
      if (alternative == "two.sided"){
        pval.asymp <- pnorm(abs(z), lower.tail = FALSE) * 2
      }else if (alternative == "less"){
        if (stat > 0){
          pval.asymp <- 1
        }else{
          pval.asymp <- pnorm(abs(z), lower.tail = FALSE)
        }
      }else if (alternative == "greater"){
        if (stat < 0){
          pval.asymp <- 1
        }else{
          pval.asymp <- pnorm(z, lower.tail = FALSE)
        }
      }
    }

    #create Monte Carlo confidence interval
    if (do.CI){
      if (!is.null(seed)){set.seed(seed)}
      beta.mc <- NA
      for (i in 1:nsims.mc){
        xy.sample <- sample(n, n, replace = TRUE)
        y.sample <- y[xy.sample]
        x.sample <- x[xy.sample] - median(x[xy.sample])
        outer.y <- outer(y.sample, y.sample, "-")
        outer.x <- outer(x.sample, x.sample, "-")
        bvals <- outer.y / outer.x
        if (do.abbreviated){
          bvals2 <- NULL
          if (n %% 2 == 0){ #even
            for (j in 1:(n / 2)){
              bvals2 <- c(bvals2, bvals[j, j + (n / 2)])
            }
          }else{ #odd
            for (j in 1:((n - 1) / 2)){
              bvals2 <- c(bvals2, bvals[j, j + ((n + 1) / 2)])
            }
          }
        }else{
          bvals2 <- bvals[upper.tri(bvals)]
        }
        #is.finite used to ignore tied x values giving infinite slopes
        beta.mc[i] <- median(bvals2[is.finite(bvals2)])
      }
      CI.mc.lower <- quantile(beta.mc, (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      CI.mc.upper <- quantile(beta.mc, 1 - (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      CI.mc.note <- "Confidence interval is basic bootstrap interval"
    }

    #create estimates of alpha to report in note
    if (do.alpha){
      d <- y - stat * x
      alpha1 <- median(d)
      Walsh.averages <- NULL
      for (i in 1:length(d)){
        for (j in i:length(d)){
          Walsh.averages <- c(Walsh.averages, (d[i] + d[j]) / 2)
        }
      }
      alpha2 <- median(Walsh.averages)
    }

    #check if message needed
    if (do.alpha){
      stat.note <- paste0("Estimate of alpha using median of d_i: ",
                          sprintf("%.5f", alpha1))
      stat.note <- paste0(stat.note, "\n")
      stat.note <- paste0(stat.note, "Hodges-Lehmann estimator of alpha: ",
                          sprintf("%.5f", alpha2))
      stat.note <- paste0(stat.note, "\n")
    }
    if ((!do.asymp && !do.exact && !do.mc) | is.null(H0)) {
      if (!is.null(stat.note)){
        stat.note <- paste0(stat.note, "\n")
      }
      stat.note <- paste0(stat.note,
                          "Neither exact, asymptotic nor Monte Carlo test requested")
    }else if (do.exact && n > max.exact.cases) {
      if (!is.null(stat.note)){
        stat.note <- paste0(stat.note, "\n")
      }
      stat.note <- paste0(stat.note,
                          "NOTE: Number of cases greater than current maximum ",
                          "allowed for exact calculations\n",
                          "required for exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ") so Monte ",
                          "Carlo p-value given")
    }
    if (do.asymp && do.abbreviated && !is.null(H0)){
      if (!is.null(stat.note)){
        stat.note <- paste0(stat.note, "\n")
      }
      stat.note <- paste0(stat.note, "NOTE: Asymptotic test not available ",
                          "when abbreviated Theil procedure requested")
    }

    #create hypotheses
    if (!is.null(H0)){
      H0val <- H0
      H0 <- paste0("H0: Theil-Kendall beta for ", varname1, " is ", H0val)
      if (alternative == "two.sided"){
        H0 <- paste0(H0, "\nH1: Theil-Kendall beta for ", varname1, " is not ",
                     H0val)
      }else if(alternative == "greater"){
        H0 <- paste0(H0, "\nH1: Theil-Kendall beta for ", varname1,
                     " is greater than ", H0val)
      }else{
        H0 <- paste0(H0, "\nH1: Theil-Kendall beta for ", varname1,
                     " is less than ", H0val)
      }
      H0 <- paste0(H0, "\n")
    }

    #return
    result <- list(title = "Theil-Kendall beta", varname1 = varname1,
                   varname2 = varname2, varname3 = varname3, stat = stat,
                   statlabel = statlabel, H0 = H0,
                   alternative = alternative, cont.corr = cont.corr, pval = pval,
                   pval.stat = pval.stat, pval.note = pval.note,
                   pval.exact = pval.exact, pval.exact.stat = pval.exact.stat,
                   pval.exact.note = pval.exact.note, targetCIwidth = CI.width,
                   actualCIwidth.exact = actualCIwidth.exact,
                   CI.exact.lower = CI.exact.lower,
                   CI.exact.upper = CI.exact.upper, CI.exact.note = CI.exact.note,
                   pval.asymp = pval.asymp, pval.asymp.stat = pval.asymp.stat,
                   pval.asymp.note = pval.asymp.note,
                   CI.asymp.lower = CI.asymp.lower,
                   CI.asymp.upper = CI.asymp.upper, CI.asymp.note = CI.asymp.note,
                   pval.mc = pval.mc, pval.mc.stat = pval.mc.stat,
                   nsims.mc = nsims.mc, pval.mc.note = pval.mc.note,
                   CI.mc.lower = CI.mc.lower, CI.mc.upper = CI.mc.upper,
                   CI.mc.note = CI.mc.note, CI.sample.lower = CI.sample.lower,
                   CI.sample.upper = CI.sample.upper, CI.sample.note = CI.sample.note,
                   stat.note = stat.note)
    class(result) <- "ANSMstat"
    return(result)
  }
