#' Perform Fisher exact test
#'
#' @description
#' `fishertest.ANSM()` is a wrapper for fisher.test() from the `stats` package - performs the Fisher exact test and is used in chapters 6, 12 and 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector or factor
#' @param y Numeric vector or factor
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10000`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.7 from "Applied Nonparametric Statistical Methods" (5th edition)
#' fishertest.ANSM(ch6$males, ch6$females)
#'
#' # Exercise 13.10 from "Applied Nonparametric Statistical Methods" (5th edition)
#' fishertest.ANSM(ch13$laid.off, ch13$employee.ages)
#'
#' @importFrom stats complete.cases median
#' @export
fishertest.ANSM <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           max.exact.cases = 10000, do.exact = TRUE) {
    stopifnot((is.vector(x) && is.numeric(x) | is.factor(x) && nlevels(x) > 1),
              (is.vector(y) && is.numeric(y) | is.factor(y) && nlevels(y) > 1),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.exact) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #default outputs
    cont.corr <- FALSE
    do.asymp <- FALSE
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
    CI.width <- NULL
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
    if (!is.factor(x) && !is.factor(y)){
      x <- x[complete.cases(x)] #remove missing cases
      y <- y[complete.cases(y)] #remove missing cases
      x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    }else{
      complete.cases.id <- complete.cases(x, y)
      x <- x[complete.cases.id] #remove missing cases
      y <- y[complete.cases.id] #remove missing cases
      if (is.factor(x)){
        x <- droplevels(x)
      }else{
        x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      }
      if (is.factor(y)){
        y <- droplevels(y)
      }else{
        y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      }
    }
    if (!is.factor(x) && !is.factor(y)){
      if (!is.null(H0)) {
        xy <- c(x - H0, y)
        varname1 <- paste0(varname1, " - ", H0)
      }else{
        H0 <- 0
        xy <- c(x, y)
      }
      med <-median(xy)
      x.gt <- sum(x > med)
      y.gt <- sum(y > med)
      nx <- x.gt + sum(x < med) #ignoring x == med
      ny <- y.gt + sum(y < med) #ignoring y == med
      n <- nx + ny
      x.mat <- matrix(c(x.gt, nx - x.gt, y.gt, ny - y.gt), nrow = 2, ncol = 2)
    }else{
      x.mat <- table(x, y)
      n <- sum(x.mat)
    }

    #Exact test
    if (do.exact && n <= max.exact.cases){
      pval.exact <- fisher.test(x.mat, alternative = alternative)$p.value
    }

    #check if message needed
    if (!do.exact) {
      test.note <- paste("Exact test not requested")
    }else if (n > max.exact.cases) {
      if (do.exact){
        test.note <- paste0("NOTE: Number of useful cases greater than ",
                            "current maximum allowed for exact\ncalculations ",
                            "required for exact test (max.exact.cases = ",
                            sprintf("%1.0f", max.exact.cases), ")")
      }
    }

    #define hypotheses
    if (!is.factor(x) && !is.factor(y)){
      if (alternative == "two.sided"){
        H0 <- paste0("H0: samples are from populations with the same median\n",
                     "H1: samples are from populations with different medians\n")
      }else if (alternative == "less"){
        H0 <- paste0("H0: samples are from populations with the same median\n",
                     "H1: median of ", varname1, " is less than median of ",
                     varname2, "\n")
      }else if (alternative == "greater"){
        H0 <- paste0("H0: samples are from populations with the same median\n",
                     "H1: median of ", varname1, " is greater than median of ",
                     varname2, "\n")
      }
    }else{
      if (alternative == "two.sided" | nlevels(x) > 2 | nlevels(y) > 2){
        H0 <- paste0("H0: ", varname1, " and ", varname2, " are independent\n",
                     "H1: ", varname1, " and ", varname2, " are not independent\n")
      }else if (alternative == "less"){
        H0 <- paste0("H0: odds ratio = 0\n",
                     "H1: odds ratio < 0\n")
      }else if (alternative == "greater"){
        H0 <- paste0("H0: odds ratio = 0\n",
                     "H1: odds ratio > 0\n")
      }
    }

    #create title
    if (!is.factor(x) && !is.factor(y)){
      title <- "Fisher exact test for median"
    }else{
      title <- "Fisher exact test"
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
