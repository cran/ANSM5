#' Perform Normal Scores test
#'
#' @description
#' `normal.scores.test()` performs the Normal Scores test and is used in chapters 6 and 8 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `25`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 5.8 from "Applied Nonparametric Statistical Methods" (5th edition)
#' normal.scores.test(ch6$groupA, ch6$groupB, do.exact = FALSE, do.asymp = TRUE)
#'
#' # Exercise 6.15 from "Applied Nonparametric Statistical Methods" (5th edition)
#' normal.scores.test(ch6$doseI, ch6$doseII)
#'
#' @importFrom stats complete.cases var
#' @importFrom utils combn
#' @export
normal.scores.test <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           max.exact.cases = 25, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

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
    y <- y[complete.cases(y)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    n.x <- length(x)
    n.y <- length(y)
    if (!is.null(H0)) {
      xy <- c(x - H0, y)
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
      xy <- c(x, y)
    }
    n.xy <- length(xy)
    #allocate ranks and calculate statistics
    xyranks <- rank(xy, ties.method = "average")
    vdW <- qnorm(xyranks / (n.xy + 1))
    vdW <- round(vdW, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    vdW.x <- sum(vdW[1:n.x])
    vdW.y <- sum(vdW[(n.x + 1):n.xy])
    if (vdW.x < vdW.y){
      n.s <- n.x
      vdW.s <- vdW.x
    }else{
      n.s <- n.y
      vdW.s <- vdW.y
    }

    #give asymptotic output if exact not possible
    if (do.exact && n.xy > max.exact.cases){
      do.asymp <- TRUE
    }

    #exact p-value
    if (do.exact && n.xy <= max.exact.cases){
      if (alternative == "two.sided"){
        pval.exact.stat <- vdW.s
        all.combn <- combn(n.xy, n.s)
        count <- 0
        for (i in 1:dim(all.combn)[2]){
          if (sum(vdW[all.combn[,i]]) <= vdW.s) {
            count <- count + 2
          }
        }
      }else if (alternative == "less"){
        pval.exact.stat <- vdW.x
        all.combn <- combn(n.xy, n.x)
        count <- 0
        for (i in 1:dim(all.combn)[2]){
          if (sum(vdW[all.combn[,i]]) <= vdW.x) {
            count <- count + 1
          }
        }
      }else if (alternative == "greater"){
        pval.exact.stat <- vdW.x
        all.combn <- combn(n.xy, n.x)
        count <- 0
        for (i in 1:dim(all.combn)[2]){
          if (sum(vdW[all.combn[,i]]) >= vdW.x) {
            count <- count + 1
          }
        }
      }
      pval.exact <- count / dim(all.combn)[2]
    }

    #asymptotic p-value (https://stat.ethz.ch/pipermail/r-help/2004-March/047190.html)
    if (do.asymp){
      if (alternative == "two.sided"){
        pval.asymp.stat <- vdW.s
        test.mean <- n.s * mean(vdW)
        test.var <- n.s * (1 - n.s / (n.xy - 1)) * (var(vdW) * (n.xy - 1) / n.xy)
        pval.asymp <- pnorm((pval.asymp.stat - test.mean) / sqrt(test.var),
                            lower.tail = TRUE) * 2
      }else{
        pval.asymp.stat <- vdW.x
        test.mean <- n.x * mean(vdW)
        test.var <- n.x * (1 - n.x / (n.xy - 1)) * (var(vdW) * (n.xy - 1) / n.xy)
        if (alternative == "greater"){
          pval.asymp <- pnorm((pval.asymp.stat - test.mean) / sqrt(test.var),
                              lower.tail = FALSE)
        }else if (alternative == "less"){
          pval.asymp <- pnorm((pval.asymp.stat - test.mean) / sqrt(test.var),
                              lower.tail = TRUE)
        }
      }
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: samples have different variances\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: variance of ", varname1, " is less than variance of ",
                   varname2, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: variance of ", varname1, " is greater than variance of ",
                   varname2, "\n")
    }

    #check if message needed
    if (do.exact && n.xy > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact\n calculations required ",
                          "(max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ")")
    }

    #return
    result <- list(title = "Normal Scores test", varname1 = varname1,
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
