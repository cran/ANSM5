#' Calculate Blomqvist coefficient
#'
#' @description
#' `blomqvist()` calculates the Blomqvist coefficient and is used in chapter 10 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector of same length as y
#' @param y Numeric vector of same length as x
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMstat object with the results from applying the function
#' @examples
#' # Example 10.9 from "Applied Nonparametric Statistical Methods" (5th edition)
#' blomqvist(ch10$q1, ch10$q2, alternative = "greater")
#'
#' # Exercise 10.7 from "Applied Nonparametric Statistical Methods" (5th edition)
#' blomqvist(ch10$ERA, ch10$SSS)
#'
#' @importFrom stats complete.cases median fisher.test
#' @export
blomqvist <-
  function(x, y, alternative = c("two.sided", "less", "greater"),
           max.exact.cases = 1000, nsims.mc = 100000, seed = NULL,
           do.exact = TRUE, do.mc = FALSE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              length(x) == length(y),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              is.logical(do.exact) == TRUE, is.logical(do.mc) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))
    varname3 <- NULL

    #unused arguments
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
    CI.sample.lower <- NULL
    CI.sample.upper <- NULL
    CI.sample.note <- NULL
    stat.note <- NULL

    #prepare
    x <- x[complete.cases(x)] #remove missing cases
    y <- y[complete.cases(y)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    n <- length(x)
    med.x <- median(x)
    med.y <- median(y)
    a <- sum(x > med.x & y > med.y)
    b <- sum(x < med.x & y > med.y)
    c <- sum(x > med.x & y < med.y)
    d <- sum(x < med.x & y < med.y)
    stat <- ((a * d) - (b * c)) / sqrt((a + b) * (c + d) * (a + c) * (b + d))
    statlabel <- "Blomqvist coefficient"

    #give mc output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.mc <- TRUE
    }

    #exact p-value
    if (do.exact && n <= max.exact.cases){
      x.mat <- matrix(c(a, b, c, d), nrow = 2, ncol = 2)
      pval.exact <- fisher.test(x.mat, alternative = alternative)$p.value
    }

    #Monte Carlo p-value
    if(do.mc){
      if (!is.null(seed)){set.seed(seed)}
      pval.mc <- 0
      for (i in 1:nsims.mc){
        x.sample <- x[sample(n, n, replace = FALSE)]
        a <- sum(x.sample > med.x & y > med.y)
        b <- sum(x.sample < med.x & y > med.y)
        c <- sum(x.sample > med.x & y < med.y)
        d <- sum(x.sample < med.x & y < med.y)
        cor.tmp <- ((a * d) - (b * c)) / sqrt((a + b) * (c + d) * (a + c) * (b + d))
        if (alternative == "two.sided"){
          if (abs(cor.tmp) >= abs(stat)){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }else if (alternative == "less"){
          if (cor.tmp <= stat){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }else if (alternative == "greater"){
          if (cor.tmp >= stat){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }
      }
    }

    #check if message needed
    if (!do.exact && !do.mc) {
      stat.note <- paste("Neither exact nor Monte Carlo test requested")
    }else if (do.exact && n > max.exact.cases) {
      stat.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact calculations\nrequired for ",
                          "exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ") so Monte ",
                          "Carlo p-value given")
    }

    #create hypotheses
    H0 <- paste0("H0: Blomqvist coefficient for ", varname1, " and ",
                 varname2, " is 0")
    if (alternative == "two.sided"){
      H0 <- paste0(H0, "\nH1: Blomqvist coefficient for ", varname1, " and ",
                   varname2, " is not 0")
    }else if(alternative == "greater"){
      H0 <- paste0(H0, "\nH1: Blomqvist coefficient for ", varname1, " and ",
                   varname2, " is greater than 0")
    }else{
      H0 <- paste0(H0, "\nH1: Blomqvist coefficient for ", varname1, " and ",
                   varname2, " is less than 0")
    }
    H0 <- paste0(H0, "\n")

    #return
    result <- list(title = "Blomqvist coefficient", varname1 = varname1,
                   varname2 = varname2, varname3 = varname3, stat = stat,
                   statlabel = statlabel, H0 = H0,
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
                   CI.mc.note = CI.mc.note, CI.sample.lower = CI.sample.lower,
                   CI.sample.upper = CI.sample.upper, CI.sample.note = CI.sample.note,
                   stat.note = stat.note)
    class(result) <- "ANSMstat"
    return(result)
  }
