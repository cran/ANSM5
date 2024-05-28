#' Perform Siegel-Tukey test
#'
#' @description
#' `siegel.tukey()` performs the Siegel-Tukey test using mean or median shift and is used in chapter 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param mean.shift Boolean indicating whether mean shift to be used instead of median shift (defaults to `FALSE`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Exercise 6.11 from "Applied Nonparametric Statistical Methods" (5th edition)
#' siegel.tukey(ch6$typeA, ch6$typeB, mean.shift = TRUE)
#'
#' # Exercise 6.16 from "Applied Nonparametric Statistical Methods" (5th edition)
#' siegel.tukey(ch6$travel, ch6$politics)
#'
#' @importFrom stats complete.cases median
#' @export
siegel.tukey <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           mean.shift = FALSE, cont.corr = TRUE, max.exact.cases = 1000,
           seed = NULL, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(seed) | is.null(seed),
              is.logical(mean.shift) == TRUE, is.logical(cont.corr) == TRUE,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #unused arguments
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
    nsims.mc <- NULL
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
    x <- x[complete.cases(x)] #remove missing cases
    y <- y[complete.cases(y)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    #equalise medians/means
    if (mean.shift){
      x <- x + (mean(y) - mean(x))
    }else{
      x <- x + (median(y) - median(x))
    }
    n.x <- length(x)
    if (!is.null(H0)) {
      xy <- c(x - H0, y)
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
      xy <- c(x, y)
    }
    n.xy <- length(xy)
    #allocate ranks
    rankstoallocate <- c(1, rep(NA, n.xy - 1))
    slbin <- 2
    slbincount <- 0
    ns <- 1
    nl <- 0
    i <- 1
    repeat{
      i <- i + 1
      if (slbincount ==2){ #swap between smallest and largest and reset count
        slbincount <- 0
        if (slbin == 1){
          slbin <- 2
        }else{
          slbin <- 1
        }
      }
      #allocate next rank
      slbincount <- slbincount + 1
      if (slbin == 1){
        ns <- ns + 1
        rankstoallocate[ns] <- i
      }else{
        nl <- nl + 1
        rankstoallocate[n.xy - nl + 1] <- i
      }
      if (i == n.xy){break}
    }
    allocatedranks <- rankstoallocate[rank(xy)]
    allocatedranksx <- allocatedranks[1:n.x]
    allocatedranksy <- allocatedranks[(n.x + 1):n.xy]

    #carry out test
    res <- wilcoxon.mann.whitney(x = allocatedranksx, y = allocatedranksy,
                                 H0 = NULL, alternative=alternative,
                                 cont.corr = cont.corr,
                                 max.exact.cases = max.exact.cases, seed = seed,
                                 do.asymp = do.asymp, do.exact = do.exact,
                                 do.CI = FALSE)
    pval <- res$pval
    pval.stat <- res$pval.stat
    pval.note <- res$pval.note
    pval.exact <- res$pval.exact
    pval.exact.stat <- res$pval.exact.stat
    pval.exact.stat <- gsub("allocatedranksx", varname1, pval.exact.stat)
    pval.exact.stat <- gsub("allocatedranksy", varname2, pval.exact.stat)
    pval.exact.note <- res$pval.exact.note
    targetCIwidth <- res$CI.width
    actualCIwidth.exact <- res$actualCIwidth.exact
    CI.exact.lower <- res$CI.exact.lower
    CI.exact.upper <- res$CI.exact.upper
    CI.exact.note <- res$CI.exact.note
    pval.asymp <- res$pval.asymp
    pval.asymp.stat <- res$pval.asymp.stat
    pval.asymp.stat <- gsub("allocatedranksx", varname1, pval.asymp.stat)
    pval.asymp.stat <- gsub("allocatedranksy", varname2, pval.asymp.stat)
    pval.asymp.note <- res$pval.asymp.note
    CI.asymp.lower <- res$CI.asymp.lower
    CI.asymp.upper <- res$CI.asymp.upper
    CI.asymp.note <- res$CI.asymp.note
    pval.mc <- res$pval.mc
    pval.mc.stat <- res$pval.mc.stat
    nsims.mc <- res$nsims.mc
    pval.mc.note <- res$pval.mc.note
    CI.mc.lower <- res$CI.mc.lower
    CI.mc.upper <- res$CI.mc.upper
    CI.mc.note <- res$CI.mc.note
    test.note <- res$test.note

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: samples have different variances\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: variance of ", varname1, " is less than variance of ",
                   varname2, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: variance of ", varname1, " is greater than variance of ",
                   varname2, "\n")
    }

    if (mean.shift){
      title <- "Siegel-Tukey test using mean shift"
    }else{
      title <- "Siegel-Tukey test using median shift"
    }

    #return
    result <- list(title = title, varname1 = varname1,
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
