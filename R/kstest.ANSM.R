#' Perform Smirnov test and Kolgomorov test
#'
#' @description
#' `kstest.ANSM()` is a wrapper for ks.test() from the `stats` package - performs the Smirnov test and Kolgomorov test and is used in chapters 4, 6 and 9 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector or a character string naming a cumulative distribution function or an actual cumulative distribution function
#' @param ... For the default method of `ks.test`, parameters of the distribution specified (as a character string) by y. Otherwise, further arguments to be passed to or from methods
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Exercise 4.3 from "Applied Nonparametric Statistical Methods" (5th edition)
#' kstest.ANSM(ch4$visiting.supporters, "pexp", rate = 2600)
#'
#' # Exercise 9.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' kstest.ANSM(ch9$boys.toothtime, ch9$girls.toothtime)
#'
#' @importFrom stats complete.cases ks.test
#' @export
kstest.ANSM <-
  function(x, y, ..., alternative = c("two.sided", "less", "greater"),
           max.exact.cases = 1000, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x),
              (is.vector(y) & is.numeric(y)) | is.character(y),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    if (is.numeric(y)){
      varname2 <- deparse(substitute(y))
    }else{
      varname2 <- shQuote(y)
    }

    #unused arguments
    cont.corr <- NULL
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
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    if (is.numeric(y)){
      y <- y[complete.cases(y)] #remove missing cases
      y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      n <- length(x) + length(y)
    }else{
      varname2 <- shQuote(y)
      n <- length(x)
    }

    #give asymptotic output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.asymp <- TRUE
    }

    #create hypotheses
    H0 <- paste0("H0: distribution of ", varname1, " matches that of ",
                 varname2)
    if (alternative == "two.sided"){
      H0 <- paste0(H0, "\nH1: distributions differ")
    }else if(alternative == "greater"){
      H0 <- paste0(H0, "\nH1: distribution of ", varname1,
                   " lies above that of ", varname2)
    }else{
      H0 <- paste0(H0, "\nH1: distribution of ", varname1,
                   " lies below that of ", varname2)
    }
    H0 <- paste0(H0, "\n")

    #exact p-value
    if(do.exact && n <= max.exact.cases){
      test.result <- ks.test(x = x, y = y, ...,
                             alternative = alternative, exact = TRUE)
      if (test.result$exact){
        pval.exact.stat <- test.result$statistic[[1]]
        pval.exact <- test.result$p.value
      }
    }

    #asymptotic p-value
    if(do.asymp){
      test.result <- ks.test(x = x, y = y, ...,
                             alternative = alternative, exact = FALSE)
      if (is.numeric(test.result$p.value)){
        pval.asymp.stat <- test.result$statistic[[1]]
        pval.asymp <- test.result$p.value
      }
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test requested")
    }else if (do.exact && n > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact calculations\nrequired for ",
                          "exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ")")
    }

    if (is.numeric(y)){
      title <- paste0("Smirnov test")
    }else{
      title <- paste0("Kolmogorov test - comparison with ", varname2)
    }

    #return
    result <- list(title = title,
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
