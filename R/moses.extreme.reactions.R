#' Perform Moses test for extreme reactions
#'
#' @description
#' `moses.extreme.reactions()` performs the Moses test for extreme reactions and is used in chapter 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.14 from "Applied Nonparametric Statistical Methods" (5th edition)
#' moses.extreme.reactions(ch6$groupI.amended, ch6$groupII)
#' moses.extreme.reactions(ch6$groupI.amended, ch6$groupII)
#'
#' @importFrom stats complete.cases
#' @export
moses.extreme.reactions <-
  function(x, y, H0 = NULL, max.exact.cases = 1000, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.exact) == TRUE)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #unused arguments
    alternative <- NULL
    cont.corr <- NULL
    CI.width <- NULL
    do.asymp <- FALSE
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
    n.x <- length(x)
    n.y <- length(y)
    n.xy <- n.x + n.y
    if (!is.null(H0)) {
      xy <- c(x - H0, y)
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
      xy <- c(x, y)
    }
    span.x <- n.xy - sum(y > max(x)) - sum(y < min(x))
    span.y <- n.xy - sum(x > max(y)) - sum(x < min(y))
    if (span.x < span.y){
      n.m <- n.x
      n.n <- n.y
      pval.exact.stat <- span.x
    }else{
      n.m <- n.y
      n.n <- n.x
      pval.exact.stat <- span.y
      tmp.obj <- varname1
      varname1 <- varname2
      varname2 <- tmp.obj
    }

    #exact p-value
    if (do.exact && n.xy <= max.exact.cases){
      pval.exact <- 0
      for (i in 0:(pval.exact.stat - n.m)){
        pval.exact <- pval.exact + choose(i + n.m - 2, i) *
          choose(n.n + 1 - i, n.n - i)
      }
      pval.exact <- pval.exact / choose(n.xy, n.m)
    }

    #define hypotheses
    H0 <- paste0("H0: extreme values are equally likely to occur in both ",
                 "populations\n",
                 "H1: extreme values are less likely to occur in ", varname1,
                 " than in ", varname2, "\n")

    #check if message needed
    if (do.exact && n.xy > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact\n calculations required ",
                          "(max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ")")
    }

    #return
    result <- list(title = "Moses test for extreme reactions",
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
