#' Calculate Pearson beta
#'
#' @description
#' `pearson.beta()` calculates the Pearson beta and is used in chapter 11 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Numeric vector of same length as x
#' @param x Numeric vector of same length as y
#' @param H0 Null hypothesis value (defaults to `NULL`)
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
#' # Example 11.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' pearson.beta(ch11$reportedtime, ch11$parentlimit, H0 = 1)
#' pearson.beta(ch11$reportedtime[1:6], ch11$parentlimit[1:6], H0 = 1)
#'
#' @importFrom stats complete.cases lm quantile
#' @export
pearson.beta <-
  function(y, x, H0 = NULL, alternative = c("two.sided", "less", "greater"),
           CI.width = 0.95, max.exact.cases = 10, nsims.mc = 100000,
           seed = NULL, do.asymp = FALSE, do.exact = TRUE, do.CI = FALSE,
           do.mc = FALSE) {
    stopifnot(is.numeric(y), is.numeric(x), length(y) == length(x),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
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
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    model <- lm(y ~ x)
    y <- model$model[, 1]
    x <- model$model[, 2]
    complete.cases.id <- complete.cases(x, y)
    y <- y[complete.cases.id] #remove missing cases
    x <- x[complete.cases.id] #remove missing cases
    n <- length(y)
    stat <- model$coefficients[2][[1]]
    statlabel <- "Pearson beta"

    if (!is.null(H0)){
      pearson.test <- pearson(x, y - H0 * x, alternative = alternative,
                              max.exact.cases = max.exact.cases,
                              nsims.mc = nsims.mc, seed = seed,
                              do.asymp = do.asymp, do.exact = do.exact,
                              do.mc = do.mc)
      pval <- pearson.test$pval
      pval.stat <- pearson.test$pval.stat
      pval.note <- pearson.test$pval.note
      pval.asymp <- pearson.test$pval.asymp
      pval.asymp.stat <- pearson.test$stat
      pval.asymp.note <- pearson.test$pval.asymp.note
      pval.exact <- pearson.test$pval.exact
      pval.exact.stat <- pearson.test$stat
      pval.exact.note <- pearson.test$pval.exact.note
      pval.mc <- pearson.test$pval.mc
      pval.mc.stat <- pearson.test$stat
      pval.mc.note <- pearson.test$pval.mc.note
      stat.note <- pearson.test$stat.note
    }

    #create Monte Carlo confidence interval
    if (do.CI){
      if (!is.null(seed)){set.seed(seed)}
      beta.mc <- NA
      for (i in 1:nsims.mc){
        xy.sample <- sample(n, n, replace = TRUE)
        y.sample <- y[xy.sample]
        x.sample <- x[xy.sample]
        beta.mc[i] <- lm(y.sample ~ x.sample)$coefficients[2]
      }
      CI.mc.lower <- quantile(beta.mc, (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      CI.mc.upper <- quantile(beta.mc, 1 - (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      CI.mc.note <- "Confidence interval is basic bootstrap interval"
    }

    #create hypotheses
    if (!is.null(H0)){
      H0val <- H0
      H0 <- paste0("H0: Pearson beta for ", varname1, " is ", H0val)
      if (alternative == "two.sided"){
        H0 <- paste0(H0, "\nH1: Pearson beta for ", varname1, " is not ", H0val)
      }else if(alternative == "greater"){
        H0 <- paste0(H0, "\nH1: Pearson beta for ", varname1,
                     " is greater than ", H0val)
      }else{
        H0 <- paste0(H0, "\nH1: Pearson beta for ", varname1, " is less than ",
                     H0val)
      }
      H0 <- paste0(H0, "\n")
    }

    #return
    result <- list(title = "Pearson beta", varname1 = varname1,
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
