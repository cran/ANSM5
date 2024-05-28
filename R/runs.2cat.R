#' Perform Runs test for two categories
#'
#' @description
#' `runs.2cat()` performs the Runs test for two categories and is used in chapters 4, 5 and 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Vector with two unique values
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.14 from "Applied Nonparametric Statistical Methods" (5th edition)
#' runs.2cat(ch4$tosses1, do.exact = FALSE, do.asymp = TRUE)
#'
#' # Exercise 6.17 from "Applied Nonparametric Statistical Methods" (5th edition)
#' runs.2cat(ch6$twins, alternative = "greater")
#'
#' @importFrom stats complete.cases pnorm
#' @export
runs.2cat <-
  function(x,  alternative=c("two.sided", "less", "greater"), cont.corr = TRUE,
           do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), length(unique(x)) == 2,
              is.logical(cont.corr) == TRUE, is.logical(do.asymp) == TRUE,
              is.logical(do.exact) == TRUE)
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

    #statistics
    x <- x[complete.cases(x)] #remove missing cases
    if (is.numeric(x)){
      x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    }
    n <- length(x)
    n1 <- sum(x == unique(x)[1])
    n2 <- sum(x == unique(x)[2])
    nruns <- 1
    for (i in 2:n){
      if (x[i] != x[i-1]){
        nruns <- nruns + 1
      }
    }

    #exact p-value
    if (do.exact){
      pval.exact.stat <- nruns
      if (alternative == "two.sided"){
        pval.exact <- NULL
        #find cumulative prob associated with nruns
        if (nruns > n / 2){
          fromlimit <- nruns
          tolimit <- n
        }else{
          fromlimit <- 1
          tolimit <- nruns
        }
        cumup1 <- 0
        for (i in fromlimit:tolimit){
          if (i %% 2 ==0){ #even
            s <- i / 2
            pr <- 2 * choose(n1 - 1, s - 1) * choose(n2 - 1, s - 1) / choose(n, n1)
            cumup1 <- cumup1 + pr
          }else{ #odd
            s <- (i - 1) / 2
            pr <- (choose(n1 - 1, s - 1) * choose(n2 - 1, s) +
                     choose(n1 - 1, s) * choose(n2 - 1, s - 1)) / choose(n, n1)
            cumup1 <- cumup1 + pr
          }
        }
        #find cumulative prob associated with opposite side
        if (nruns > n / 2){
          fromlimit <- 1
          tolimit <- nruns - 1
        }else{
          fromlimit <- n
          tolimit <- nruns + 1
        }
        cumup2 <- 0
        for (i in fromlimit:tolimit){
          if (i %% 2 ==0){ #even
            s <- i / 2
            pr <- 2 * choose(n1 - 1, s - 1) * choose(n2 - 1, s - 1) / choose(n, n1)
            cumup2 <- cumup2 + pr
          }else{ #odd
            s <- (i - 1) / 2
            pr <- (choose(n1 - 1, s - 1) * choose(n2 - 1, s) +
                     choose(n1 - 1, s) * choose(n2 - 1, s - 1)) / choose(n, n1)
            cumup2 <- cumup2 + pr
          }
          if (cumup2 - cumup1 > 1e-07){
            pval.exact <- cumup1 + cumup2 - pr
            break
          }
        }
        if (is.null(pval.exact)){pval.exact <- cumup1 + cumup2} #will sum to 1
      }else if (alternative == "less"){
        pval.exact <- 0
        for (i in 1:nruns){
          if (i %% 2 ==0){ #even
            s <- i / 2
            pr <- 2 * choose(n1 - 1, s - 1) * choose(n2 - 1, s - 1) / choose(n, n1)
            pval.exact <- pval.exact + pr
          }else{ #odd
            s <- (i - 1) / 2
            pr <- (choose(n1 - 1, s - 1) * choose(n2 - 1, s) +
                     choose(n1 - 1, s) * choose(n2 - 1, s - 1)) / choose(n, n1)
            pval.exact <- pval.exact + pr
          }
        }
      }else if (alternative == "greater"){
        pval.exact <- 0
        for (i in nruns:n){
          if (i %% 2 ==0){ #even
            s <- i / 2
            pr <- 2 * choose(n1 - 1, s - 1) * choose(n2 - 1, s - 1) / choose(n, n1)
            pval.exact <- pval.exact + pr
          }else{ #odd
            s <- (i - 1) / 2
            pr <- (choose(n1 - 1, s - 1) * choose(n2 - 1, s) +
                     choose(n1 - 1, s) * choose(n2 - 1, s - 1)) / choose(n, n1)
            pval.exact <- pval.exact + pr
          }
        }
      }
    }

    #asymptotic p-value
    if (do.asymp){
      ExpR <- 1 + 2 * n1 * n2 / n
      VarR <- (2 * n1 * n2 * (2 * n1 * n2 - n)) / (n * n * (n - 1))
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
      if (n < 20){
        pval.asymp.note <-
          pval.asymp.note <- paste("(WARNING: n is less than 20 so asymptotic",
                                   "test is not recommended)")
      }
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test requested")
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
    result <- list(title = "Runs test for two categories" , varname1 = varname1,
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
