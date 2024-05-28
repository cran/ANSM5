#' Perform Log odds ratio test
#'
#' @description
#' `logoddsratio.2x2()` performs the Log odds ratio test and is used in chapter 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Binary factor of same length as y
#' @param y Binary factor of same length as x
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Exercise 13.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' #logoddsratio.2x2(ch13$physical.activity[ch13$gender == "Boy"],
#' #  ch13$tv.viewing[ch13$gender == "Boy"], do.exact = FALSE, do.asymp = TRUE)
#' #logoddsratio.2x2(ch13$physical.activity[ch13$gender == "Girl"],
#' #  ch13$tv.viewing[ch13$gender == "Girl"], do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases r2dtable pnorm
#' @export
logoddsratio.2x2 <-
  function(x, y, max.exact.cases = 10, nsims.mc = 100000,
           seed = NULL, do.exact = TRUE, do.asymp = FALSE, do.mc = FALSE) {
    stopifnot(is.factor(x), is.factor(y), nlevels(x) == 2, nlevels(y) == 2,
              length(x) == length(y), all(table(x, y) > 0),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.exact) == TRUE, is.logical(do.asymp) == TRUE,
              is.logical(do.mc) == TRUE)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #unused arguments
    H0 <- NULL
    cont.corr <- NULL
    alternative <- NULL
    CI.width <- NULL
    do.CI <- FALSE
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
    test.note <- NULL

    #prepare
    complete.cases.id <- complete.cases(x, y)
    x <- x[complete.cases.id] #remove missing cases
    y <- y[complete.cases.id] #remove missing cases
    x <- droplevels(x)
    y <- droplevels(y)
    n <- length(x)
    tab.n <- nlevels(x) * nlevels(y)
    rtots <- table(x)
    ctots <- table(y)
    tab.xy <- table(x, y)
    stat <- log((tab.xy[1, 1] * tab.xy[2, 2]) / (tab.xy[1, 2] * tab.xy[2, 1]),
                base = exp(1))

    #give mc output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.mc <- TRUE
    }

    #exact p-value
    if(do.exact && n <= max.exact.cases){
      pval.exact.stat <- stat
      permutations <- perms(n)
      n.perms <- dim(permutations)[1]
      pval.exact <- 0
      for (i in 1:n.perms){
        tmp.tab.xy <- table(x[permutations[i,]], y)
        tmp.stat <- log((tmp.tab.xy[1, 1] * tmp.tab.xy[2, 2]) /
                          (tmp.tab.xy[1, 2] * tmp.tab.xy[2, 1]), base = exp(1))
        if (tmp.stat >= pval.exact.stat){
          pval.exact <- pval.exact + 2 / n.perms
        }
      }
      if (pval.exact > 1){pval.exact <- 1}
    }

    #Monte Carlo p-value
    if (do.mc){
      pval.mc.stat <- stat
      if (!is.null(seed)){set.seed(seed)}
      pval.mc <- 0
      for (i in 1:nsims.mc){
        tmp.tab.xy <- r2dtable(1, rtots, ctots)[[1]]
        tmp.stat <- log((tmp.tab.xy[1, 1] * tmp.tab.xy[2, 2]) /
                          (tmp.tab.xy[1, 2] * tmp.tab.xy[2, 1]), base = exp(1))
        if (tmp.stat >= pval.mc.stat){
          pval.mc <- pval.mc + 2 / nsims.mc
        }
      }
      if (pval.mc > 1){pval.mc <- 1}
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- stat
      var <- sum(1 / tab.xy)
      pval.asymp <- pnorm(stat / sqrt(var), lower.tail = FALSE) * 2
    }

    #check if message needed
    if (!do.exact && !do.mc && !do.asymp) {
      test.note <- paste("Neither exact, asymptotic nor Monte Carlo test requested")
    }else if (do.exact && n > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact calculations\nrequired for ",
                          "exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ") so Monte ",
                          "Carlo p-value given")
    }

    #define hypotheses
    H0 <- paste0("H0: ", varname1, " and ", varname2, " are independent\n",
                 "H1: ", varname1, " and ", varname2, " are not independent\n")

    #return
    result <- list(title = "Log odds ratio test", varname1 = varname1,
                   varname2 = varname2, H0 = H0,
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
                   CI.mc.note = CI.mc.note,
                   test.note = test.note)
    class(result) <- "ANSMtest"
    return(result)
  }
