#' Performs Lilliefors test of Normality
#'
#' @description
#' `lilliefors()` performs Lilliefors test of Normality and is used in chapters 4, 5 and 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param alternative Type of alternative hypothesis (defaults to `c("two.sided")`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.4 from "Applied Nonparametric Statistical Methods" (5th edition)
#' lilliefors(ch4$ages, seed = 1)
#'
#' # Exercise 6.15 from "Applied Nonparametric Statistical Methods" (5th edition)
#' lilliefors(ch6$doseI.2, seed = 1, nsims = 1000)
#'
#' @importFrom stats complete.cases rnorm sd
#' @export
lilliefors <-
  function(x, alternative = c("two.sided"), nsims.mc = 10000, seed = NULL) {
    stopifnot(is.vector(x), is.numeric(x),
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed))
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

    #statistics
    x <- x[complete.cases(x)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    x <- sort(x)
    n <- length(x)
    zx <- scale(x)
    zx.pnorm <- pnorm(zx)
    zx.pnorm <- unique(zx.pnorm)
    s <- rank(x, ties.method = "max") / n
    s <- unique(s)
    diff1 <- zx.pnorm - c(0, s[1:(length(s) - 1)])
    diff2 <- s - zx.pnorm

    #Monte Carlo p-value
    if (!is.null(seed)){set.seed(seed)}
    diffs.sim <- NULL
    diffs1.sim <- NULL
    diffs2.sim <- NULL
    for (i in 1:nsims.mc){
      x.sim <- rnorm(n, mean(x), sd(x))
      x.sim <- sort(x.sim)
      zx.sim <- scale(x.sim)
      zx.sim.pnorm <- pnorm(zx.sim)
      s.sim <- rank(x.sim, ties.method = "max") / n
      diff1.sim <- zx.sim.pnorm - c(0, s.sim[1:(length(s.sim) - 1)])
      diff2.sim <- s.sim - zx.sim.pnorm
      diffs.sim <- c(diffs.sim, max(abs(diff1.sim), abs(diff2.sim)))
      diffs1.sim <- c(diffs1.sim, max(abs(diff1.sim)))
      diffs2.sim <- c(diffs2.sim, max(abs(diff2.sim)))
    }
    pval.mc.stat <- max(abs(diff1), abs(diff2))
    pval.mc <- sum(pval.mc.stat < diffs.sim) / nsims.mc

    #create hypotheses
    H0 <- paste0("H0: distribution of ", varname1, " is Normal\n",
                 "H1: distribution of ", varname1, " is not Normal\n")

    #return
    result <- list(title = "Lilliefors test of Normality",
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
