#' Perform Binomial test
#'
#' @description
#' `binom()` performs the Binomial test and calculates the Binomial confidence interval and is used in chapters 4, 5 and 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param r Number of successes
#' @param n Number of trials
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10000000`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.6 from "Applied Nonparametric Statistical Methods" (5th edition)
#' binom(3, 20)
#'
#' # Exercise 5.8 from "Applied Nonparametric Statistical Methods" (5th edition)
#' binom(24, 40, 0.5)
#'
#' @importFrom stats dbinom pbinom pnorm qnorm
#' @export
binom <-
  function(r, n, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           CI.width = 0.95, max.exact.cases = 10000000, do.asymp = FALSE,
           do.exact = TRUE, do.CI = TRUE) {
    stopifnot(is.numeric(r), length(r) == 1,
              is.numeric(n), length(n) == 1, r < n,
              is.numeric(H0) | is.null(H0), length(H0) == 1 | is.null(H0),
              H0 >= 0 | is.null(H0), H0 <= 1 | is.null(H0),
              CI.width > 0, CI.width < 1,
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- paste0("r = ", r, ", n = ", n)

    #default outputs
    varname2 <- NULL
    cont.corr <- NULL
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

    #give asymptotic output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.asymp <- TRUE
    }

    #exact p-value
    if (!is.null(H0) && do.exact && n <= max.exact.cases){
      if (alternative=="two.sided"){
        if (r / n == H0){
          pval.exact <- 1
        }else if (r / n < H0){
          db <- dbinom(r, n, H0)
          for (i in n:0){
            if (dbinom(i, n, H0) > db){
              break
            }
          }
          pval.exact <- pbinom(r, n, H0) + pbinom(i, n, H0, lower.tail = FALSE)
        }else{
          db <- dbinom(r, n, H0)
          for (i in 0:n){
            if (dbinom(i, n, H0) > db){
              break
            }
          }
          pval.exact <- pbinom(i - 1, n, H0) + pbinom(r - 1, n, H0, lower.tail = FALSE)
        }
      }else if (alternative == "less"){
        pval.exact <- pbinom(r, n, H0)
      }else if (alternative == "greater"){
        pval.exact <- pbinom(r - 1, n, H0, lower.tail = FALSE)
      }
    }

    #exact CI
    if (do.exact && do.CI){
      if(alternative == "two.sided"){
        divisor <- 2
      }else{
        divisor <- 1
      }
      if (alternative == "less"){
        CI.exact.lower <- 0
      }else{
        CI.exact.lower <- NULL
        p <- 1
        repeat{
          if (1 - pbinom(r - 1, n, p) < (1 - CI.width) / divisor) {
            CI.exact.lower <- p
            break
          }
          p <- p - 0.00001
        }
      }
      if (alternative == "greater"){
        CI.exact.upper <- 1
      }else{
        CI.exact.upper <- NULL
        p <- 0
        repeat{
          if (pbinom(r, n, p) < (1 - CI.width) / divisor) {
            CI.exact.upper <- p
            break
          }
          p <- p + 0.00001
        }
      }
      if (alternative == "two.sided"){
        actualCIwidth.exact <- 1 - pbinom(r, n, CI.exact.upper) - (1 - pbinom(r - 1, n, CI.exact.lower))
      }else if (alternative == "less"){
        actualCIwidth.exact <- 1 - pbinom(r, n, CI.exact.upper)
      }else{
        actualCIwidth.exact <- pbinom(r - 1, n, CI.exact.lower)
      }
      if (is.null(CI.exact.lower) | is.null(CI.exact.upper) | is.null(actualCIwidth.exact)){
        actualCIwidth.exact <- NULL
        CI.exact.lower <- NULL
        CI.exact.upper <- NULL
      }
    }

    #asymptotic p-value
    if (!is.null(H0) && do.asymp){
      phat <- r / n
      pval.asymp.stat <- (phat - H0) / sqrt(phat * (1 - phat) / n)
      if (alternative == "two.sided"){
        if (pval.asymp.stat == 0){
          pval.asymp <- 1
        }else if (pval.asymp.stat < 0){
          pval.asymp <- pnorm(pval.asymp.stat) + pnorm(-pval.asymp.stat, lower.tail = FALSE)
        }else{
          pval.asymp <- pnorm(pval.asymp.stat, lower.tail = FALSE) + pnorm(-pval.asymp.stat)
        }
      }else if (alternative == "less"){
        if (pval.asymp.stat >= 0){
          pval.asymp <- 1
        }else{
          pval.asymp <- pnorm(pval.asymp.stat)
        }
      }else{
        if (pval.asymp.stat <= 0){
          pval.asymp <- 1
        }else{
          pval.asymp <- pnorm(pval.asymp.stat, lower.tail = FALSE)
        }
      }
      if (n < 20){
        pval.asymp.note <-
          pval.asymp.note <- paste("(WARNING: n is less than 20 so asymptotic",
                                   "test is not recommended)")
      }
    }

    #asymptotic CI
    if (do.asymp && do.CI){
      if(alternative == "two.sided"){
        divisor <- 2
      }else{
        divisor <- 1
      }
      phat <- r / n
      z1 <- qnorm((1 - CI.width) / divisor)
      z2 <- qnorm(1 - (1 - CI.width) / divisor)
      if (alternative == "less"){
        CI.asymp.lower <- 0
      }else{
        CI.asymp.lower <- phat + z1 * sqrt(phat * (1 - phat) / n)
      }
      if (alternative == "greater"){
        CI.asymp.upper <- 1
      }else{
        CI.asymp.upper <- phat + z2 * sqrt(phat * (1 - phat) / n)
      }
      if (n < 20){
        CI.asymp.note <- paste("(WARNING: n is less than 20 so asymptotic",
                               "CI is not recommended)")
      }
    }

    #check if message needed
    if ((!do.asymp && !do.exact) | !do.CI) {
      test.note <- paste("Neither exact nor asymptotic test/confidence ",
                         "interval requested")
    }else if (do.exact && n > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact calculations\nrequired for ",
                          "exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ")")
    }

    #define title
    if (is.null(H0)){
      title <- "Binomial confidence interval"
    }else{
      title <- "Binomial test"
    }

    #return
    result <- list(title = title, varname1 = varname1, varname2 = varname2,
                   H0 = H0, alternative = alternative, cont.corr = cont.corr,
                   pval = pval, pval.stat = pval.stat, pval.note = pval.note,
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
