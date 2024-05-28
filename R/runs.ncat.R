#' Perform Runs test for three or more categories
#'
#' @description
#' `runs.ncat()` performs the Runs test for three or more categories and is used in chapters 4 and 7 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Vector or factor
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `TRUE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.15 from "Applied Nonparametric Statistical Methods" (5th edition)
#' runs.ncat(ch4$births, alternative = "less")
#'
#' # Exercise 7.16 from "Applied Nonparametric Statistical Methods" (5th edition)
#' runs.ncat(ch7$regions[order(ch7$affordability)], alternative = "less")
#'
#' @importFrom stats complete.cases pnorm
#' @export
runs.ncat <-
  function(x,  alternative=c("two.sided", "less", "greater"), cont.corr = TRUE,
           nsims.mc = 100000, seed = NULL, do.asymp = TRUE, do.mc = FALSE) {
    stopifnot(is.vector(x) | is.factor(x), is.logical(cont.corr) == TRUE,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.mc) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))

    #default outputs
    varname2 <- NULL
    CI.width <- NULL
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

    #statistics
    x <- x[complete.cases(x)] #remove missing cases
    if (is.factor(x)){
      x <- as.vector(x)
    }else{
      if (is.numeric(x)){
        x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      }
    }
    n <- length(x)
    k <- length(unique(x))
    pri <- rep(NA, k)
    for (i in 1:k){
      pri[i] <- sum(x == unique(x)[i]) / n
    }
    nruns <- 1
    for (i in 2:n){
      if (x[i] != x[i-1]){
        nruns <- nruns + 1
      }
    }

    #asymptotic p-value
    if (do.asymp){
      ExpR <- n * (1 - sum(pri ^ 2)) + 1
      VarR <- n * (sum(pri ^ 2 - 2 * pri ^ 3) + sum(pri ^ 2) ^ 2)
      if (cont.corr && nruns != ExpR){
        pval.asymp.stat <- (nruns - ExpR + (0.5 - (nruns > ExpR))) / sqrt(VarR)
      }else{
        pval.asymp.stat <- (nruns - ExpR) / sqrt(VarR)
      }
      if (alternative == "two.sided"){
        if (pval.asymp.stat < 0){
          pval.asymp <- 2 * pnorm(pval.asymp.stat)
        }else{
          pval.asymp <- 2 * pnorm(-pval.asymp.stat)
        }
      }else if (alternative == "less"){
        pval.asymp <- pnorm(pval.asymp.stat)
      }else if (alternative == "greater"){
        pval.asymp <- 1 - pnorm(pval.asymp.stat)
      }
      if (n < 12){
        pval.asymp.note <-
          pval.asymp.note <- paste("(WARNING: n is less than 12 so asymptotic",
                                   "test is not recommended)")
      }
    }

    #Monte Carlo p-value
    if (do.mc){
      if (!is.null(seed)){set.seed(seed)}
      nruns.sim <- NULL
      for (i in 1:nsims.mc){
        x.sim <- sample(x)
        nr.sim <- 1
        for (j in 2:n){
          if (x.sim[j] != x.sim[j-1]){
            nr.sim <- nr.sim + 1
          }
        }
        if (nr.sim == 1){message(i, x.sim)}
        nruns.sim <- c(nruns.sim, nr.sim)
      }
      if (alternative == "two.sided"){
        pval.mc <- 2 * min(sum(nruns >= nruns.sim) / nsims.mc,
                           sum(nruns <= nruns.sim) / nsims.mc)
      }else if (alternative == "less"){
        pval.mc <- sum(nruns >= nruns.sim) / nsims.mc
      }else if (alternative == "greater"){
        pval.mc <- sum(nruns <= nruns.sim) / nsims.mc
      }
    }

    #check if message needed
    if (!do.asymp && !do.mc) {
      test.note <- paste("Neither exact nor Monte Carlo test requested")
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: number of runs consistent with randomness\n",
                   "H1: number of runs not consistent with randomness\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: number of runs consistent with randomness\n",
                   "H1: number of runs fewer than consistent with randomness\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: number of runs consistent with randomness\n",
                   "H1: number of runs greater than consistent with randomness\n")
    }

    #return
    result <- list(title = "Runs test for three or more categories" ,
                   varname1 = varname1, varname2 = varname2,
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
