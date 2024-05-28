#' Perform Hodges-Ajne test
#'
#' @description
#' `hodges.ajne()` performs the Hodges-Ajne test and is used in chapter 4 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param alternative Type of alternative hypothesis (defaults to `c("two.sided")`)
#' @param minx Minimum value for x (defaults to `0`)
#' @param maxx Maximum value for x (defaults to `360`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 4.16 from "Applied Nonparametric Statistical Methods" (5th edition)
#' hodges.ajne(ch4$times.as.degrees)
#'
#' # Exercise 4.14 from "Applied Nonparametric Statistical Methods" (5th edition)
#' hodges.ajne(ch4$board.angles)
#'
#' @importFrom stats complete.cases
#' @export
hodges.ajne <-
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

    #calculate angles to check
    angles <- NULL
    check3 <- 0
    can_exit <- FALSE
    starta <- x[1]
    repeat{
      if (check3 == 0){
        if (length(x[x > starta]) > 0){
          nexta <- min(x[x > starta])
        }else{
          nexta <- x[1]
        }
      }else{
        if (length(x[x > check1]) > 0){
          nexta <- min(x[x > check1]) - 180 + 360 * (min(x[x > check1]) - 180 < 0)
        }else{
          nexta <- x[1] - 180 + 360 * (x[1] - 180 < 0)
        }
      }
      check1 <- starta + 180 - 360 * (starta + 180 > 360)
      check2 <- nexta + 180 - 360 * (nexta + 180 > 360)
      check3 <- sum(x > check1 & x < check2)
      if (check3 == 0){
        if (starta < nexta){
          angles <- c(angles, (starta + nexta) / 2)
        }else{
          tempa <- ((starta - 360) + nexta) / 2
          if (tempa < 0){
            angles <- c(angles, 360 + tempa)
          }else{
            angles <- c(angles, tempa)
          }
        }
        starta <- nexta
        can_exit <- TRUE
      }
      if (starta == x[1] && can_exit){break}
    }

    #calculate m
    pval.stat <- length(x) + 1
    for (i in 1:length(angles)){
      m.a1 <- sum(x - angles[i] + 360 * (x - angles[i] < 0) > 180)
      m.a2 <- sum(x - angles[i] + 360 * (x - angles[i] < 0) < 180)
      if (m.a1 < pval.stat){pval.stat <- m.a1}
      if (m.a2 < pval.stat){pval.stat <- m.a2}
    }
    #cope with all angles in one half of circle
    for (i in 1:length(x)){
      if (i < length(x)){
        if (x[i + 1] - x[i] >= 180){
          pval.stat <- 0
          break
        }
      }else{
        if (x[i] - x[1] <= 180){
          pval.stat <- 0
          break
        }
      }
    }

    #p-value
    if (pval.stat < n / 3){
      pval <- choose(n, pval.stat) * (n - 2 * pval.stat) / (2 ^ (n - 1))
    }else{
      pval <- NA
      pval.note <- paste0(varname1, " is distributed in a sufficiently\n",
                          "uniform fashion to result in the smallest number\n",
                          "of cases found in a half circle being more than a\n",
                          "third of the total number of cases, meaning that\n",
                          "the Mardia (1972) formula to calculate the\n",
                          "probabilities cannot be used")
    }

    #create hypotheses
    H0 <- paste0("H0: distribution of ", varname1, " is uniform\n",
                 "H1: distribution of ", varname1, " is not uniform\n")

    #return
    result <- list(title = "Hodges-Ajne test",
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
