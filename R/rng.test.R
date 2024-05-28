#' Perform Range test
#'
#' @description
#' `rng.test()` performs the Range test and is used in chapter 4 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param alternative Type of alternative hypothesis (defaults to `c("two.sided")`)
#' @param minx Minimum value for x (defaults to `0`)
#' @param maxx Maximum value for x (defaults to `360`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.17 from "Applied Nonparametric Statistical Methods" (5th edition)
#' rng.test(ch4$dates.as.degrees)
#'
#' # Exercise 4.13 from "Applied Nonparametric Statistical Methods" (5th edition)
#' rng.test(ch4$accident.bearings)
#'
#' @importFrom stats complete.cases
#' @export
rng.test <-
  function(x, alternative = c("two.sided"), minx = 0, maxx = 360) {
    stopifnot(is.vector(x), is.numeric(x), length(x) > 1)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))

    #default outputs
    varname2 <- NULL
    cont.corr <- NULL
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

    #prepare
    x <- x[complete.cases(x)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    n <- length(x)
    x <- (x + minx) #make smallest possible number zero
    x <- x * (360 / maxx) #make largest possible number 360
    x <- sort(x)

    #calculate largest arc containing all cases
    arclength <- 360 - x[length(x)] + x[1] #for length overlapping 360
    for (i in 2:length(x)){
      if (x[i] - x[i - 1] > arclength){
        arclength <- x[i] - x[i - 1]
      }
    }
    pval.stat <- (360 - arclength) * pi / 180

    #p-value
    pval <- 0
    for (k in 1:n){
      kpart <- 1 - (k * (2 * pi - pval.stat) / (2 * pi))
      if (kpart < 0){break}
      pval <- pval + ((-1) ^ (k - 1)) * choose(n, k) * (kpart ^ (n - 1))
    }

    #create hypotheses
    H0 <- paste0("H0: distribution of ", varname1, " is uniform\n",
                 "H1: distribution of ", varname1, " is not uniform\n")

    #return
    result <- list(title = "Range test",
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
