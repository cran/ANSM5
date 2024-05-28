#' Perform Cramer-von Mises test
#'
#' @description
#' `cramer.von.mises()` performs the Cramer-von Mises test and is used in chapter 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.16 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cramer.von.mises(ch6$salivaF, ch6$salivaM)
#' cramer.von.mises(ch6$salivaF, ch6$salivaM, alternative = "greater")
#'
#' @importFrom stats complete.cases ecdf
#' @export
cramer.von.mises <-
  function(x, y, alternative = c("two.sided", "less", "greater")) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y))
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #unused arguments
    max.exact.cases <- NULL
    cont.corr <- NULL
    CI.width <- NULL
    do.exact <- FALSE
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
    cdf_x <- ecdf(x)(x)
    cdf_y <- ecdf(y)(y)
    cdf_xy <- c(cdf_x, cdf_y)[order(rank(c(x, y)))]
    xy_labs <- c(rep("x", length(x)), rep("y", length(y)))[order(rank(c(x, y)))]
    if (xy_labs[1] == "x"){
      Sx <- cdf_xy[1]
      Sy <- 0
    }else{
      Sx <- 0
      Sy <- cdf_xy[1]
    }
    for (i in 2:length(xy_labs)){
      if (xy_labs[i] == "x"){
        Sx <- c(Sx, cdf_xy[i])
        Sy <- c(Sy, Sy[i - 1])
      }else{
        Sx <- c(Sx, Sx[i - 1])
        Sy <- c(Sy, cdf_xy[i])
      }
    }

    #create hypotheses
    H0 <- paste0("H0: distribution of ", varname1, " matches that of ",
                 varname2)
    if (alternative == "two.sided"){
      H0 <- paste0(H0, "\nH1: distributions differ")
    }else if(alternative == "greater"){
      H0 <- paste0(H0, "\nH1: distribution of ", varname1,
                   " lies above that of ", varname2)
    }else{
      H0 <- paste0(H0, "\nH1: distribution of ", varname1,
                   " lies below that of ", varname2)
    }
    H0 <- paste0(H0, "\n")

    #p-value
    S2d <- sum((Sx - Sy) ** 2)
    pval.stat <- length(x) * length (y) * S2d / ((length(x) + length(y)) ** 2)
    if (alternative == "two.sided"){
      if (pval.stat > 0.743){
        pval <- "p-value < 0.01"
      }else if (pval.stat > 0.461){
        pval <- "p-value < 0.05"
      }else if (pval.stat > 0.347){
        pval <- "p-value < 0.10"
      }else{
        pval <- "p-value >= 0.10"
      }
    }else if (alternative == "greater"){
      if (sum(Sx > Sy) < sum(Sx < Sy)){ #"less" more likely than "greater"
        pval <- paste0("Greater evidence for opposite of alternative ",
        "hypothesis than for alternative hypothesis")
      }else{
        if (pval.stat > 0.743){
          pval <- "p-value < 0.005"
        }else if (pval.stat > 0.461){
          pval <- "p-value < 0.025"
        }else if (pval.stat > 0.347){
          pval <- "p-value < 0.05"
        }else{
          pval <- "p-value >= 0.05"
        }
      }
    }else if (alternative == "less"){
      if (sum(Sx > Sy) > sum(Sx < Sy)){ #"greater" more likely than "less"
        pval <- paste0("Greater evidence for opposite of alternative ",
                       "hypothesis than for alternative hypothesis")
      }else{
        if (pval.stat > 0.743){
          pval <- "p-value < 0.005"
        }else if (pval.stat > 0.461){
          pval <- "p-value < 0.025"
        }else if (pval.stat > 0.347){
          pval <- "p-value < 0.05"
        }else{
          pval <- "p-value >= 0.05"
        }
      }
    }

    #return
    result <- list(title = "Cramer-von Mises test",
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
