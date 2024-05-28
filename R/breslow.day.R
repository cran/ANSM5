#' Perform Breslow and Day test
#'
#' @description
#' `breslow.day()` performs the Breslow and Day test and is used in chapter 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Binary factor of same length as y, z
#' @param y Binary factor of same length as x, z
#' @param z Factor of same length as x, y
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 13.3 from "Applied Nonparametric Statistical Methods" (5th edition)
#' breslow.day(ch13$machine, ch13$output.status, ch13$material.source)
#'
#' # Exercise 13.7 from "Applied Nonparametric Statistical Methods" (5th edition)
#' breslow.day(ch13$medicine, ch13$response, ch13$location)
#'
#' @importFrom stats complete.cases pchisq qnorm
#' @export
breslow.day <-
  function(x, y, z, CI.width = 0.95, do.asymp = TRUE, do.CI = TRUE) {
    stopifnot(is.factor(x), is.factor(y), is.factor(z),
              nlevels(x) == 2, nlevels(y) == 2, nlevels(z) > 1,
              length(x) == length(y), length(y) == length(z),
              length(CI.width) == 1, CI.width > 0, CI.width < 1,
              is.logical(do.asymp) == TRUE, is.logical(do.CI) == TRUE)

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
    #calculate theta
    theta.numer <- 0
    theta.denom <- 0
    tab <- array(rep(NA, nlevels(z) * 2 * 2), dim = c(nlevels(z), 2, 2))
    n1plus.z <- rep(NA, nlevels(z))
    nplus1.z <- rep(NA, nlevels(z))
    N.z <- rep(NA, nlevels(z))
    for (i in 1:nlevels(z)){
      tab[i,,] <- table(x[z == levels(z)[i]], y[z == levels(z)[i]])
      theta.numer <- theta.numer + tab[i, 1, 1] * tab[i, 2, 2] / sum(tab[i,,])
      theta.denom <- theta.denom + tab[i, 1, 2] * tab[i, 2, 1] / sum(tab[i,,])
      n1plus.z[i] <- rowSums(tab[i,,])[1]
      nplus1.z[i] <- colSums(tab[i,,])[1]
      N.z[i] <- sum(tab[i,,])
    }
    theta <- theta.numer / theta.denom
    #calculate Breslow and Day M2
    M.BD <- 0
    for (i in 1:nlevels(z)){
      a = 1 - theta
      b = N.z[i] - n1plus.z[i] - nplus1.z[i] + (theta * (n1plus.z[i] + nplus1.z[i]))
      c = -1 * theta * n1plus.z[i] * nplus1.z[i]
      x.1 <- (-b + (sqrt(b ^ 2 - 4 * a * c))) / (2 * a)
      x.2 <- (-b - (sqrt(b ^ 2 - 4 * a * c))) / (2 * a)
      if (x.1 > n1plus.z[i] | x.1 > nplus1.z[i]){
        m11 <- x.2
      }else{
        m11 <- x.1
      }
      m12 <- rowSums(tab[i,,])[1] - m11
      m21 <- colSums(tab[i,,])[1] - m11
      m22 <- sum(tab[i,,]) - m11 - m12 - m21
      M.BD <- M.BD + (m11 - tab[i, 1, 1]) ^ 2 /
        (1 / (1 / m11 + 1 / m12 + 1 / m21 + 1 / m22))
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- M.BD
      pval.asymp <- pchisq(pval.asymp.stat, nlevels(z) - 1, lower.tail = FALSE)
      #confidence intervals
      if (do.CI){
        a <- rep(NA, nlevels(z))
        b <- rep(NA, nlevels(z))
        c <- rep(NA, nlevels(z))
        d <- rep(NA, nlevels(z))
        ac <- rep(NA, nlevels(z))
        ad <- rep(NA, nlevels(z))
        bc <- rep(NA, nlevels(z))
        bd <- rep(NA, nlevels(z))
        for (i in 1:nlevels(z)){
          a[i] <- (tab[i, 1, 1] + tab[i, 2, 2]) / sum(tab[i,,])
          b[i] <- (tab[i, 1, 2] + tab[i, 2, 1]) / sum(tab[i,,])
          c[i] <- (tab[i, 1, 1] * tab[i, 2, 2]) / sum(tab[i,,])
          d[i] <- (tab[i, 1, 2] * tab[i, 2, 1]) / sum(tab[i,,])
          ac[i] <- a[i] * c[i]
          ad[i] <- a[i] * d[i]
          bc[i] <- b[i] * c[i]
          bd[i] <- b[i] * d[i]
        }
        P <- sum(ac) / (2 * sum(c) ^ 2)
        Q <- sum(ad + bc) / (2 * sum(c) * sum(d))
        R <- sum(bd) / (2 * sum(d) ^ 2)
        Var <- P + Q + R
        CI.asymp.lower <- exp(log(theta) -
          qnorm((1 - CI.width) / 2, lower.tail = FALSE) * sqrt(Var))
        CI.asymp.upper <- exp(log(theta) +
          qnorm((1 - CI.width) / 2, lower.tail = FALSE) * sqrt(Var))
      }
    }

    #check if message needed
    if (!do.mc && !do.asymp) {
      test.note <- paste("Asymptotic test not requested")
    }

    #define hypotheses
    H0 <- paste0("H0: Odds ratio for ", varname1, " by ", varname2,
                 " is the same for all levels of ", varname3, "\n",
                 "H1: Odds ratio for ", varname1, " by ", varname2,
                 " is not the same for all levels of ", varname3, "\n")

    #return
    result <- list(title = "Breslow and Day test",
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
