#' Perform Hettmansperger and Elmore interaction test
#'
#' @description
#' `hettmansperger.elmore()` performs the Hettmansperger and Elmore interaction test and is used in chapter 8 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Numeric vector of same length as factor.a, factor.b
#' @param factor.a Factor of same length as y, factor.b
#' @param factor.b Factor of same length as y, factor.a
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `1000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `TRUE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @param median.polish Boolean indicating whether or not to use median polish (defaults to `FALSE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 8.6 from "Applied Nonparametric Statistical Methods" (5th edition)
#' hettmansperger.elmore(ch8$plant.weight.2, ch8$growth.hormone, ch8$undersoil.heating)
#'
#' # Exercise 8.3 from "Applied Nonparametric Statistical Methods" (5th edition)
#' hettmansperger.elmore(ch8$game.time, ch8$experience, ch8$game)
#'
#' @importFrom stats complete.cases aov anova medpolish resid pchisq
#' @export
hettmansperger.elmore <-
  function(y, factor.a, factor.b, nsims.mc = 1000, seed = NULL,
           do.asymp = TRUE, do.mc = FALSE, median.polish = FALSE) {
    stopifnot(is.vector(y), is.numeric(y), is.factor(factor.a), is.factor(factor.b),
              length(y) == length(factor.a), length(factor.a) == length(factor.b),
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.mc) == TRUE,
              is.logical(median.polish) == TRUE)

    #labels
    varname1 <- deparse(substitute(y))
    varname2 <- paste0(deparse(substitute(factor.a)), ", ",
                       deparse(substitute(factor.b)))

    #unused arguments
    alternative <- NULL
    cont.corr <- NULL
    CI.width <- NULL
    do.CI <- FALSE
    do.exact <- FALSE
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
    complete.cases.ID <- complete.cases(y, factor.a, factor.b)
    y <- y[complete.cases.ID] #remove missing cases
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    factor.a <- factor.a[complete.cases.ID] #remove missing cases
    factor.b <- factor.b[complete.cases.ID] #remove missing cases
    nlev.a <- nlevels(factor.a)
    nlev.b <- nlevels(factor.b)

    #calculate statistics
    if (median.polish){
      tmp.y <- y
      repeat{
        medians.a <- tapply(tmp.y, factor.a, median)
        tmp.y <- tmp.y - medians.a[as.numeric(factor.a)]
        medians.b <- tapply(tmp.y, factor.b, median)
        tmp.y <- tmp.y - medians.b[as.numeric(factor.b)]
        if (sum(medians.a, medians.b) < 0.0001){break}
      }
      resids.to.rank <- tmp.y
    }else{
      y.aov <- aov(y ~ factor.a + factor.b)
      resids.to.rank <- resid(y.aov)
    }
    tmp.anova <- anova(aov(rank(resids.to.rank) ~ factor.a + factor.b
                           + factor.a * factor.b))
    int.ss <- tmp.anova$`Sum Sq`[3]
    tot.ms <- sum(tmp.anova$`Sum Sq`) / sum(tmp.anova$Df)
    Q <- int.ss / tot.ms

    #asymptotic p-value
    if(do.asymp){
      pval.asymp.stat <- Q
      pval.asymp <- pchisq(pval.asymp.stat, (nlev.a - 1) * (nlev.b - 1),
                           lower.tail = FALSE)
    }

    #Monte Carlo p-value
    if (do.mc){
      pval.mc.stat <- Q
      if (!is.null(seed)){set.seed(seed)}
      pval.mc <- 0
      for (i in 1:nsims.mc){
        y_i <- sample(y, length(y))
        if (median.polish){
          tmp.y <- y_i
          repeat{
            medians.a <- tapply(tmp.y, factor.a, median)
            tmp.y <- tmp.y - medians.a[as.numeric(factor.a)]
            medians.b <- tapply(tmp.y, factor.b, median)
            tmp.y <- tmp.y - medians.b[as.numeric(factor.b)]
            if (sum(medians.a, medians.b) < 0.0001){break}
          }
          resids.to.rank_i <- tmp.y
        }else{
          y.aov_i <- aov(y_i ~ factor.a + factor.b)
          resids.to.rank_i <- resid(y.aov_i)
        }
        y.aov_i <- aov(y_i ~ factor.a + factor.b)
        tmp.anova_i <- anova(aov(rank(resids.to.rank_i) ~ factor.a + factor.b
                               + factor.a * factor.b))
        int.ss <- tmp.anova_i$`Sum Sq`[3]
        tot.ms <- sum(tmp.anova_i$`Sum Sq`) / sum(tmp.anova_i$Df)
        stat_i <- int.ss / tot.ms
        if (stat_i >= pval.mc.stat){
          pval.mc <- pval.mc + 1 / nsims.mc
        }
      }
    }

    #check if message needed
    if (!do.asymp && !do.mc) {
      test.note <- paste("Neither asymptotic nor Monte Carlo test requested")
    }

    #create hypotheses
    H0 <- paste0("H0: there is no interaction between ", varname2, "\n")
    H0 <- paste0(H0, "H1: there is an interaction between ", varname2)
    H0 <- paste0(H0, "\n")

    #return
    result <- list(title = "Hettmansperger and Elmore interaction test",
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
