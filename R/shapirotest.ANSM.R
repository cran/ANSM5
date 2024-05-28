#' Perform Shapiro-Wilk test of Normality
#'
#' @description
#' `shapirotest.ANSM()` is a wrapper for shapiro.test() from the `stats` package - performs the Shapiro-Wilk test of Normality and is used in chapters 4 and 5 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param alternative Type of alternative hypothesis (defaults to `c("two.sided")`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.4 from "Applied Nonparametric Statistical Methods" (5th edition)
#' shapirotest.ANSM(ch4$ages)
#'
#' # Example 5.3 from "Applied Nonparametric Statistical Methods" (5th edition)
#' shapirotest.ANSM(ch5$bp.incorrect)
#'
#' @importFrom stats complete.cases shapiro.test
#' @export
shapirotest.ANSM <-
  function(x,  alternative = c("two.sided")) {
    stopifnot(is.vector(x), is.numeric(x))
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))

    #default outputs
    varname2 <- NULL
    cont.corr <- NULL
    CI.width <- NULL
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
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues

    #create hypotheses
    H0 <- paste0("H0: distribution of ", varname1, " is Normal\n",
                 "H1: distribution of ", varname1, " is not Normal\n")

    #p-value
    test.result <- shapiro.test(x = x)
    if (is.numeric(test.result$p.value)){
      pval.asymp.stat <- test.result$statistic[[1]]
      pval.asymp <- test.result$p.value
    }

    #return
    result <- list(title = "Shapiro-Wilk test of Normality",
                   varname1 = varname1, varname2 = varname2, H0 = H0,
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
