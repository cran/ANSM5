#' Perform Gehan-Wilcoxon test
#'
#' @description
#' `gehan.wilcoxon()` performs the Gehan-Wilcoxon test and is used in chapter 9 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector of same length as y, x.c, y.c
#' @param y Numeric vector of same length as x, x.c, y.c
#' @param x.c Binary vector of same length as x, y, x.c
#' @param y.c Binary vector of same length as x, y, y.c
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `100000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 9.1 from "Applied Nonparametric Statistical Methods" (5th edition)
#' gehan.wilcoxon(ch9$symp.survtime, ch9$asymp.survtime,
#'   ch9$symp.censor, ch9$asymp.censor, alternative = "less",
#'   do.exact = FALSE, do.asymp = TRUE)
#'
#' # Exercise 9.5 from "Applied Nonparametric Statistical Methods" (5th edition)
#' gehan.wilcoxon(ch9$regimeA.survtime, ch9$regimeB.survtime,
#'   ch9$regimeA.censor, ch9$regimeB.censor, do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases pnorm
#' @importFrom utils combn
#' @export
gehan.wilcoxon <-
  function(x, y, x.c, y.c, alternative=c("two.sided", "less", "greater"),
           max.exact.perms = 100000, nsims.mc = 100000, seed = NULL,
           do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              is.vector(x.c), is.numeric(x.c), is.vector(y.c), is.numeric(y.c),
              all(x.c == 0 | x.c == 1), all(y.c == 0 | y.c == 1),
              length(x) == length(x.c), length(y) == length(y.c),
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
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
    x <- x[complete.cases(x, x.c)] #remove missing cases
    x.c <- x.c[complete.cases(x, x.c)] #remove missing cases
    y <- y[complete.cases(y, y.c)] #remove missing cases
    y.c <- y.c[complete.cases(y, y.c)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    xy <- c(x, y)
    xy.c <- c(x.c, y.c)
    n <- length(xy)
    n.perms <- choose(n, min(length(x), length(y)))
    #calculate statistics
    U.x <- 0
    for (i in 1:length(x)){
      if (x.c[i] == 0){
        U.x <- U.x + sum(y > x[i]) + sum((y == x[i]) * y.c) +
          0.5 * sum(y == x[i] * (1 - y.c)) + 0.5 * sum((y < x[i]) * y.c)
      }else{
        U.x <- U.x + sum((y > x[i]) * (1 - y.c)) * 0.5 + 0.5 * sum(y.c)
      }
    }
    U.y <- 0
    for (i in 1:length(y)){
      if (y.c[i] == 0){
        U.y <- U.y + sum(x > y[i]) + sum((x == y[i]) * x.c) +
          0.5 * sum(x == y[i] * (1 - x.c)) + 0.5 * sum((x < y[i]) * x.c)
      }else{
        U.y <- U.y + sum((x > y[i]) * (1 - x.c)) * 0.5 + 0.5 * sum(x.c)
      }
    }
    if (alternative == "two.sided"){
      U <- max(U.x, U.y)
    }else if (alternative == "less"){
      U <- U.x
    }else{
      U <- U.y
    }

    #exact p-value
    if (do.exact){
      OverflowState <- FALSE
      if (n.perms <= max.exact.perms){
        #try complete combinations
        try_result <- suppressWarnings(try(
          combins <- combn(n, min(length(x), length(y))), silent = TRUE)
        )
        if (any(class(try_result) == "try-error")){
          OverflowState <- TRUE
        }else{
          n.combins <- dim(combins)[2]
        }
      }
      if (n.perms > max.exact.perms | OverflowState){
        #use Monte Carlo
        if (!is.null(seed)){set.seed(seed)}
        n.combins <- nsims.mc
      }
      #evaluate all combinations
      tmp.pval <- 0
      for (i in 1:n.combins){
        #define data
        if (n.perms > max.exact.perms | OverflowState){
          tmp.combins <- sample(n, min(length(x), length(y)))
        }else{
          tmp.combins <- combins[,i]
        }
        tmp.x <- xy[tmp.combins]
        tmp.x.c <- xy.c[tmp.combins]
        tmp.y <- xy[-tmp.combins]
        tmp.y.c <- xy.c[-tmp.combins]
        #calculate test statistics
        tmp.U.x <- 0
        for (i in 1:length(tmp.x)){
          if (tmp.x.c[i] == 0){
            tmp.U.x <- tmp.U.x + sum(tmp.y > tmp.x[i]) +
              sum((tmp.y == tmp.x[i]) * tmp.y.c) +
              0.5 * sum(tmp.y == tmp.x[i] * (1 - tmp.y.c)) +
              0.5 * sum((tmp.y < tmp.x[i]) * tmp.y.c)
          }else{
            tmp.U.x <- tmp.U.x + sum((tmp.y > tmp.x[i]) * (1 - tmp.y.c)) * 0.5 +
              0.5 * sum(tmp.y.c)
          }
        }
        tmp.U.y <- 0
        for (i in 1:length(tmp.y)){
          if (tmp.y.c[i] == 0){
            tmp.U.y <- tmp.U.y + sum(tmp.x > tmp.y[i]) +
              sum((tmp.x == tmp.y[i]) * tmp.x.c) +
              0.5 * sum(tmp.x == tmp.y[i] * (1 - tmp.x.c)) +
              0.5 * sum((tmp.x < tmp.y[i]) * tmp.x.c)
          }else{
            tmp.U.y <- tmp.U.y + sum((tmp.x > tmp.y[i]) * (1 - tmp.x.c)) * 0.5 +
              0.5 * sum(tmp.x.c)
          }
        }
        if (alternative == "two.sided"){
          tmp.U <- max(tmp.U.x, tmp.U.y)
        }else if (alternative == "less"){
          tmp.U <- tmp.U.x
        }else{
          tmp.U <- tmp.U.y
        }
        #check against actual test statistic
        if (tmp.U > U){
          tmp.pval <- tmp.pval + 1 / n.combins
        }
      }
      #output
      if (n.perms <= max.exact.perms && !OverflowState){
        pval.exact.stat <- U
        pval.exact <- tmp.pval
      }else{
        pval.mc.stat <- U
        pval.mc <- tmp.pval
      }
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- U
      n.x <- length(x)
      n.y <- length(y)
      z <- (pval.asymp.stat - (n.x * n.y / 2)) /
        sqrt(n.x * n.y * (n.x + n.y + 1) / 12)
      pval.asymp <- pnorm(z, lower.tail = FALSE)
      if (alternative == "two.sided"){
        pval.asymp <- pval.asymp * 2
      }
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test requested")
    }else if (do.exact && n.perms > max.exact.perms) {
      test.note <- paste0("NOTE: Number of permutations required greater than ",
                          "current maximum allowed for exact calculations\n",
                          "required for exact test (max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ") so Monte ",
                          "Carlo p-value given")
    }else if (do.exact && OverflowState){
      test.note <- paste0("NOTE: Insufficient memory for exact calculations",
                          "required for exact test\n(max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ") so Monte ",
                          "Carlo p-value given")
    }
    if (do.asymp){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      test.note <- paste("NOTE: asymptotic test not adjusted for ties")
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: the survival times distributions are identical\n",
                   "H1: the survival times distribution of ", varname1,
                   " are different from those of ", varname2, "\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: the survival times distributions are identical\n",
                   "H1: the survival times distribution of ", varname1,
                   " are below those of ", varname2, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: the survival times distributions are identical\n",
                   "H1: the survival times distribution of ", varname1,
                   " are above those of ", varname2, "\n")
    }

    #return
    result <- list(title = "Gehan-Wilcoxon test", varname1 = varname1,
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
