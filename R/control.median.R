#' Perform Control median test
#'
#' @description
#' `control.median()` performs the Control median test and is used in chapters 6 and 9 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.9 from "Applied Nonparametric Statistical Methods" (5th edition)
#' control.median(ch6$sampleI, ch6$sampleII, alternative = "greater")
#'
#' # Exercise 9.8 from "Applied Nonparametric Statistical Methods" (5th edition)
#' control.median(ch9$bulbA, ch9$bulbB, alternative = "greater", nsims = 1000)
#'
#' @importFrom stats complete.cases median pnorm
#' @export
control.median <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           CI.width = 0.95, max.exact.cases = 1000,
           nsims.mc = 10000, seed = NULL,
           do.asymp = FALSE, do.exact = TRUE, do.CI = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              CI.width > 0, CI.width < 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- paste0(deparse(substitute(x)), " (Control)")
    varname2 <- deparse(substitute(y))

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
    test.note <- NULL

    #prepare
    x <- x[complete.cases(x)] #remove missing cases
    y <- y[complete.cases(y)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    nxy <- length(x) + length(y)
    if (!is.null(H0)) {
      x <- x - H0
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
    }
    med_x <- median(x)

    #give asymptotic output if exact not possible
    if (do.exact && nxy > max.exact.cases){
      do.asymp <- TRUE
    }

    #exact p-value
    if (do.exact && nxy <= max.exact.cases){
      #test statistic
      if (alternative == "two.sided"){
        pval.exact.stat <- min(sum(y < med_x), sum(y > med_x))
      }else if (alternative == "greater"){
        pval.exact.stat <- sum(y < med_x)
      }else{
        pval.exact.stat <- sum(y > med_x)
      }
      #calculate
      m <- length(x)
      n <- length(y)
      if (m %% 2 == 0){
        r <- m / 2
      }else{
        r <- (m + 1) / 2
      }
      pval.exact <- 0
      if (alternative == "two.sided" | alternative == "greater"){
        for (i in 0:pval.exact.stat){
          pval.exact <- pval.exact +
            choose(r + i, i) * choose(m - r + n - i, n - i) /
            choose(m + n + 1, n)
        }
      }
      if (alternative == "two.sided" | alternative == "less"){
        for (i in seq(n, n - pval.exact.stat, -1)){
          pval.exact <- pval.exact +
            choose(r + i, i) * choose(m - r + n - i, n - i) /
            choose(m + n + 1, n)
        }
      }
    }

    #asymptotic p-value
    if (do.asymp){
      #test statistic
      if (alternative == "two.sided"){
        pval.asymp.stat <- min(sum(y < med_x), sum(y > med_x))
      }else if (alternative == "greater"){
        pval.asymp.stat <- sum(y < med_x)
      }else{
        pval.asymp.stat <- sum(y > med_x)
      }
      #calculate (N.B. m and n switched from exact test)
      m <- length(y)
      n <- length(x)
      Z <- (pval.asymp.stat - (m / 2)) / sqrt(m * (m + n) / (4 * n))
      pval.asymp <- pnorm(Z, lower.tail = TRUE)
      if (alternative == "two.sided"){pval.asymp <- pval.asymp * 2}
    }

    #confidence interval
    if (do.CI){
      bs.ci.res <- bs(x = x, y = y, CI.width = CI.width, nsims.bs = nsims.mc,
                         seed = seed)$CI
      CI.mc.lower <- bs.ci.res[1]
      CI.mc.upper <- bs.ci.res[2]
      CI.mc.note <- paste0("Confidence interval for difference (", varname1,
      " minus ", varname2, ")\nis basic bootstrap interval for the median")
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test requested")
    }else if (do.exact && nxy > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact\ncalculations required ",
                          "for exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ")")
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: samples are from populations with the same median\n",
                   "H1: samples are from populations with different medians\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: samples are from populations with the same median\n",
                   "H1: median of ", varname2, " is less than median of ",
                   varname1, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: samples are from populations with the same median\n",
                   "H1: median of ", varname2, " is greater than median of ",
                   varname1, "\n")
    }

    #return
    result <- list(title = "Control median test", varname1 = varname1,
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
