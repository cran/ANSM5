#' Perform Sign test
#'
#' @description
#' `sgn.test()` performs the Sign test and is used in chapters 3, 4, 5 and 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector, or binary factor and H0 is NULL
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000000`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 3.1 from "Applied Nonparametric Statistical Methods" (5th edition)
#' sgn.test(ch3$sampleI, 110)
#'
#' # Exercise 6.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' sgn.test(ch5$LVF - ch5$RVF, 0)
#'
#' @importFrom stats complete.cases dbinom pbinom pnorm qnorm
#' @export
sgn.test <-
  function(x, H0 = NULL, alternative = c("two.sided", "less", "greater"),
           cont.corr = TRUE, CI.width = 0.95, max.exact.cases = 1000000,
           do.asymp = FALSE, do.exact = TRUE, do.CI = TRUE) {
  stopifnot((is.vector(x) && is.numeric(x) | is.factor(x)),
            if(is.factor(x)){is.null(H0) && nlevels(x) == 2}else{TRUE},
            ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
            is.numeric(max.exact.cases), length(max.exact.cases) == 1,
            is.logical(cont.corr) == TRUE, is.numeric(CI.width),
            length(CI.width) == 1, CI.width > 0, CI.width < 1,
            is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
            is.logical(do.CI) == TRUE)
  alternative <- match.arg(alternative)

  #labels
  varname1 <- deparse(substitute(x))

  #default outputs
  varname2 <- NULL
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

  #statistics
  x <- x[complete.cases(x)] #remove missing cases
  if (is.factor(x)){
    pluses <- sum(as.numeric(x) > 1.5)
    H0 <- paste0("H0: there are no differences between the levels of ",
                 varname1, "\n",
                 "H1: there is a difference between the levels of ", varname1,
                 "\n")
  }else{
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    if (!is.null(H0)) {
      x <- x[x != H0] #remove cases equal to H0
      pluses <- sum(x > H0)
    }
  }
  n <- length(x)

  #give asymptotic output if exact not possible
  if (do.exact && n > max.exact.cases){
    do.asymp <- TRUE
  }

  #exact p-value
  if (!is.null(H0) && do.exact && n <= max.exact.cases){
    pval.exact.less <- pbinom(pluses, n, 0.5)
    pval.exact.greater <- pbinom(n-pluses, n, 0.5)
    if (alternative=="two.sided"){
      pval.exact <- min(1, min(pval.exact.less, pval.exact.greater) * 2)
    }else if (alternative == "less"){
      pval.exact <- pval.exact.less
    }else if (alternative == "greater"){
      pval.exact <- pval.exact.greater
    }
  }

  #exact CI
  if (do.exact && do.CI && n <= max.exact.cases){
    lower <- 0
    lower.cumu <- 0
    repeat{
      if (lower.cumu + dbinom(lower, n, 0.5) > (1 - CI.width) / 2) {break}
      lower.cumu <- lower.cumu + dbinom(lower, n, 0.5)
      lower <- lower + 1
    }
    if (lower == 0){
      CI.exact.lower = Inf
    }else{
      CI.exact.lower <- sort(x, decreasing = FALSE)[lower]
    }
    upper <- n
    upper.cumu <- 0
    repeat{
      if (upper.cumu + dbinom(upper, n, 0.5) > (1 - CI.width) / 2) {break}
      upper.cumu <- upper.cumu + dbinom(upper, n, 0.5)
      upper <- upper - 1
    }
    if (upper == n){
      CI.exact.upper = Inf
    }else{
      CI.exact.upper <- sort(x, decreasing = FALSE)[upper+1]
    }
    actualCIwidth.exact <- 1 - lower.cumu - upper.cumu
    if (is.null(CI.exact.lower) | is.null(CI.exact.upper) | is.null(actualCIwidth.exact)){
      actualCIwidth.exact <- NULL
      CI.exact.lower <- NULL
      CI.exact.upper <- NULL
    }
  }

  #asymptotic p-value (with/without continuity correction)
  if (!is.null(H0) && do.asymp){
    if (alternative=="two.sided"){
      pval.asymp <- pnorm((pluses - n * 0.5 + (0.5 - (pluses > n * 0.5)) *
                             (cont.corr == TRUE))/(0.5 * sqrt(n)),
                          lower.tail = pluses < n * 0.5) * 2
    }else if (alternative == "less"){
      pval.asymp <- pnorm((pluses - n * 0.5 + (0.5 - (pluses > n * 0.5)) *
                             (cont.corr == TRUE))/(0.5 * sqrt(n)),
                          lower.tail = TRUE)
    }else if (alternative == "greater"){
      pval.asymp <- pnorm((pluses - n * 0.5 + (0.5 - (pluses > n * 0.5)) *
                             (cont.corr == TRUE))/(0.5 * sqrt(n)),
                          lower.tail = FALSE)
    }
    if (pluses < 5 | n - pluses < 5){
      pval.asymp.note <- paste("(WARNING: n * p and/or n * (1 - p) are less",
                               "than 5 so asymptotic test is not recommended)")
    }
  }

  #asymptotic CI
  if (do.asymp && do.CI){
    meanrank <- n / 2
    sdrank <- sqrt(n) / 2
    S <- qnorm((1 - CI.width) / 2) * sdrank + meanrank
    CI.asymp.lower <- sort(x, decreasing = FALSE)[round(S, 0) + 1]
    CI.asymp.upper <- sort(x, decreasing = FALSE)[n - round(S, 0)]
    if (n < 20){
      CI.asymp.note <- paste("(WARNING: n is less than 20 so asymptotic",
                             "CI is not recommended)")
    }
  }

  #check if message needed
  if (!do.asymp && !do.exact) {
    test.note <- paste("Neither exact nor asymptotic test/confidence interval ",
                       "requested")
  }else if ((do.CI | do.exact) && n > max.exact.cases) {
    test.note <- paste0("NOTE: Number of useful cases greater than current ",
                        "maximum allowed for exact calculations\nrequired for ",
                        "exact test or confidence interval (max.exact.cases = ",
                        sprintf("%1.0f", max.exact.cases), ")")
  }
  if (is.null(H0) && (do.exact | do.asymp)){
    if (!is.null(test.note)){
      test.note <- paste0(test.note, "\n")
    }
    test.note <- paste0(test.note, "NOTE: No value for H0 given so no p-value ",
                        "can be calculated")
  }


  #return
  result <- list(title = "Sign test", varname1 = varname1, varname2 = varname2,
                 H0 = H0,
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
