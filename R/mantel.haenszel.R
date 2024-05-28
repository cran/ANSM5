#' Perform Mantel-Haenszel test
#'
#' @description
#' `mantel.haenszel()` performs the Mantel-Haenszel test and is used in chapter 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Binary factor of same length as y, z
#' @param y Binary factor of same length as x, z
#' @param z Factor of same length as x, y
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 13.4 from "Applied Nonparametric Statistical Methods" (5th edition)
#' mantel.haenszel(ch13$drug, ch13$side.effects, ch13$age.group)
#'
#' #  from "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @importFrom stats complete.cases pchisq
#' @export
mantel.haenszel <-
  function(x, y, z, do.asymp = TRUE) {
    stopifnot(is.factor(x), is.factor(y), is.factor(z),
              nlevels(x) == 2, nlevels(y) == 2, nlevels(z) > 1,
              length(x) == length(y), length(y) == length(z),
              is.logical(do.asymp) == TRUE)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))
    varname3 <- deparse(substitute(z))

    #unused arguments
    H0 <- NULL
    cont.corr <- NULL
    alternative <- NULL
    do.exact <- FALSE
    do.mc <- FALSE
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
    nsims.mc <- NULL
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
    complete.cases.id <- complete.cases(x, y, z)
    x <- x[complete.cases.id] #remove missing cases
    y <- y[complete.cases.id] #remove missing cases
    z <- z[complete.cases.id] #remove missing cases
    x <- droplevels(x)
    y <- droplevels(y)
    z <- droplevels(z)
    #calculate M2 and variance
    tab <- array(rep(NA, nlevels(z) * 2 * 2), dim = c(nlevels(z), 2, 2))
    n1plus.z <- rep(NA, nlevels(z))
    nplus1.z <- rep(NA, nlevels(z))
    n2plus.z <- rep(NA, nlevels(z))
    nplus2.z <- rep(NA, nlevels(z))
    N.z <- rep(NA, nlevels(z))
    for (i in 1:nlevels(z)){
      tab[i,,] <- table(x[z == levels(z)[i]], y[z == levels(z)[i]])
      n1plus.z[i] <- rowSums(tab[i,,])[1]
      nplus1.z[i] <- colSums(tab[i,,])[1]
      n2plus.z[i] <- rowSums(tab[i,,])[2]
      nplus2.z[i] <- colSums(tab[i,,])[2]
      N.z[i] <- sum(tab[i,,])
    }
    m11.plus <- 0
    varn11.z <- rep(NA, nlevels(z))
    for (i in 1:nlevels(z)){
      m11.plus <- m11.plus + n1plus.z[i] * nplus1.z[i] / N.z[i]
      varn11.z[i] <- (n1plus.z[i] * nplus1.z[i] * n2plus.z[i] * nplus2.z[i]) /
        ((N.z[i] ^ 2) * (N.z[i] - 1))
    }
    M2.numer <- (m11.plus - sum(tab[, 1, 1])) ^ 2
    M2 <- M2.numer / sum(varn11.z)

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- M2
      pval.asymp <- pchisq(pval.asymp.stat, 1, lower.tail = FALSE)
    }

    #check if message needed
    if (!do.mc && !do.asymp) {
      test.note <- paste("Asymptotic test not requested")
    }

    #define hypotheses
    H0 <- paste0("H0: Odds ratio for ", varname1, " by ", varname2,
                 " is equal to 1 for all levels of ", varname3, "\n",
                 "H1: Odds ratio for ", varname1, " by ", varname2,
                 " is not equal to 1 for all levels of ", varname3, "\n")

    #return
    result <- list(title = "Mantel-Haenszel test",
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
