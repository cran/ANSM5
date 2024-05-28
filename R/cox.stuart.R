#' Perform Cox-Stuart test
#'
#' @description
#' `cox.stuart()` performs the Cox-Stuart test and is used in chapters 4 and 10 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10000000`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.13 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cox.stuart(ch4$precipitation)
#'
#' # Exercise 10.5 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cox.stuart(ch10$age.at.death, alternative = "less")
#'
#' @importFrom stats complete.cases
#' @export
cox.stuart <-
  function(x, alternative = c("two.sided", "less", "greater"), cont.corr = TRUE,
           max.exact.cases = 10000000, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.logical(cont.corr) == TRUE,
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)
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
    CI.width <- NULL
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
    if (length(x) %% 2 == 1){ #odd so remove middle value
      x <- c(x[1:floor(length(x)/2)], x[(ceiling(length(x)/2) +  1):length(x)])
    }
    diffs <- x[1:(length(x)/2)] - x[(length(x)/2 + 1):length(x)]

    #give asymptotic output if exact not possible
    if (do.exact && length(x) > max.exact.cases){
      do.asymp <- TRUE
    }

    #carry out test
    H0 <- 0
    do.CI <- FALSE
    result <- sgn.test(diffs, H0 = H0, alternative = alternative,
                       cont.corr = cont.corr,
                       max.exact.cases = max.exact.cases, do.asymp = do.asymp,
                       do.exact = do.exact, do.CI = do.CI)
    H0 <- paste0("H0: no monotonic trend\n",
                 "H1: montonic trend exists\n")
    pval.exact <- result$pval.exact
    pval.asymp <- result$pval.asymp
    pval.asymp.note <- result$pval.asymp.note
    test.note <- result$test.note

    #return
    result <- list(title = "Cox-Stuart test", varname1 = varname1,
                   varname2 = varname2,
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
