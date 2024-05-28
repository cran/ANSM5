#' Perform Kruskal-Wallis test with van der Waerden scores
#'
#' @description
#' `kruskal.wallis.vdW()` performs the Kruskal-Wallis test with van der Waerden scores and is used in chapter 7 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector of same length as g
#' @param g Factor of same length as x
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `15`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 7.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' kruskal.wallis.vdW(ch7$age, ch7$positions)
#' kruskal.wallis.vdW(ch7$age, ch7$positions, do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases
#' @importFrom utils combn
#' @export
kruskal.wallis.vdW <-
  function(x, g, max.exact.cases = 15, nsims.mc = 10000, seed = NULL,
           do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.factor(g), length(x) == length(g),
              length(x[complete.cases(x)]) == length(g[complete.cases(g)]),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(g))

    #unused arguments
    alternative <- NULL
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
    x <- x[complete.cases(x)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    g <- g[complete.cases(g)] #remove missing cases
    n <- length(x)
    rank.x <- rank(x)
    vdW <- qnorm(rank.x / (n + 1))
    table_g <- table(g)
    c <- 0
    sr <- sum(vdW ** 2)

    #calculate test statistic
    T <- ((n - 1) * (sum(by(vdW, g, sum) ** 2 / table_g) - c)) / (sr - c)
    T <- round(T, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues

    #exact p-value
    if (do.exact && n <= max.exact.cases){
      combins <- NULL
      for (ig in 1:(nlevels(g) - 1)){
        if (ig == 1){
          combins <- t(combn(n, table_g[ig]))
        }else{
          combins2 <- NULL
          for (i in 1:dim(combins)[1]){
            combins2 <- rbind(combins2,
                              cbind(matrix(rep(combins[i,],
                                               choose(n - dim(combins)[2],
                                                      table_g[ig])),
                                           ncol = dim(combins)[2],
                                           byrow = TRUE),
                                    t(combn(setdiff(seq(1:n), combins[i,]),
                                            table_g[ig]))))
          }
          combins <- combins2
        }
      }
      combins <- data.frame(matrix(vdW[combins], ncol=dim(combins)[2]))
      combins$T <- 0
      for (i in 1:nlevels(g)){
        if (i < nlevels(g)){
          combins$T <- combins$T +
            rowSums(combins[which(g == levels(g)[i])]) ** 2 / table_g[i]
        }else{
          combins$T <- combins$T +
            (sum(vdW) - rowSums(combins[,1:(n - table_g[i])])) ** 2 / table_g[i]
        }
      }
      combins$T <- ((n - 1) * (combins$T - c)) / (sr - c)
	  combins$T <- round(combins$T, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      pval.exact <- sum(combins$T >= T) / dim(combins)[1]
      pval.exact.stat <- T
    }

    #Monte Carlo p-value
    if (do.exact && n > max.exact.cases){
      if (!is.null(seed)){set.seed(seed)}
      T.sim <- NULL
      for (i in 1:nsims.mc){
        rank.sim <- sample(vdW, n, replace = FALSE)
        Ti <- ((n - 1) * (sum(by(rank.sim, g, sum) ** 2 / table_g) - c)) /
          (sr - c)
        T.sim <- c(T.sim, Ti)
      }
      pval.mc <- sum(T.sim >= T) / nsims.mc
      pval.mc.stat <- T
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp <- pchisq(T, nlevels(g) - 1, lower.tail = FALSE)
      pval.asymp.stat <- T
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test requested")
    }else if (n > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact\ncalculations required ",
                          "for exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ") so Monte ",
                          "Carlo p-value given")
    }

    #define hypotheses
    H0 <- paste0("H0: samples are from the same population\n",
                 "H1: samples differ in location\n")

    #return
    result <- list(title = "Kruskal-Wallis test with van der Waerden scores",
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
