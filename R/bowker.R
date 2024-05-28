#' Perform Bowker's extension of McNemar's test
#'
#' @description
#' `bowker()` performs the Bowker's extension of McNemar's test and is used in chapter 12 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Factor of same length as y, or two-dimensional square table
#' @param y Factor of same length as x (or NULL if x is table) (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 12.12 from "Applied Nonparametric Statistical Methods" (5th edition)
#' bowker(ch12$side.effect.new, ch12$side.effect.old)
#'
#' # Exercise 12.12 from "Applied Nonparametric Statistical Methods" (5th edition)
#' bowker(ch12$first.response, ch12$second.response)
#'
#' @importFrom stats complete.cases pchisq
#' @export
bowker <-
  function(x, y = NULL, do.asymp = TRUE) {
    stopifnot((is.factor(x) && nlevels(x) > 1 |
                 is.table(x) && dim(x)[1] == dim(x)[2]),
              ((is.factor(y) && nlevels(x) == nlevels(y)) |
                 (is.table(x) && is.null(y))),
              (nlevels(y) > 1 | (is.table(x) && is.null(y))),
              (length(x) == length(y) | is.null(y)),
              is.logical(do.asymp) == TRUE)

    #labels
    if (is.table(x)){
      varname1 <- names(attributes(x)$dimnames)[1]
      varname2 <- names(attributes(x)$dimnames)[2]
    }else{
      varname1 <- deparse(substitute(x))
      varname2 <- deparse(substitute(y))
    }

    #unused arguments
    H0 <- NULL
    alternative <- NULL
    cont.corr <- NULL
    nsims.mc <- NULL
    do.exact <- NULL
    do.mc <- NULL
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
    if (is.table(x)){
      y <- NULL
      n <- sum(x)
      tab <- x
    }else{
      complete.cases.id <- complete.cases(x, y)
      x <- x[complete.cases.id] #remove missing cases
      x <- droplevels(x)
      y <- y[complete.cases.id] #remove missing cases
      y <- droplevels(y)
      n <- length(x)
      tab <- table(x, y)
    }
    n.levels <- dim(tab)[1]
    stat <- 0
    for (i in 1:(n.levels - 1)){
      for (j in (i + 1):n.levels){
        stat <- stat + ((tab[i, j] - tab[j, i]) ^ 2) / (tab[i, j] + tab[j, i])
      }
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- stat
      pval.asymp <- pchisq(stat, n.levels * (n.levels - 1) / 2,
                           lower.tail = FALSE)
    }

    #check if message needed
    if (!do.asymp) {
      test.note <- paste("No test requested")
    }

    #define hypotheses
    H0 <- paste0("H0: ", varname1, " and ", varname2, " are independent\n",
                 "H1: ", varname1, " and ", varname2, " are not independent\n")

    #return
    result <- list(title = "Bowker's extension of McNemar's test",
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
