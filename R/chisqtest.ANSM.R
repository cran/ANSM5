#' Perform Chi-squared test
#'
#' @description
#' `chisqtest.ANSM()` is a wrapper for chisq.test() from the `stats` package - performs the Chi-squared test and is used in chapters 12 and 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Factor of same length as y, or table
#' @param y Factor of same length as x (or NULL if x is table) (defaults to `NULL`)
#' @param p Vector of probabilities (expressed as numbers between 0 and 1 and summing to 1) of same length as x or NULL (defaults to `NULL`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 12.1 from "Applied Nonparametric Statistical Methods" (5th edition)
#' chisqtest.ANSM(ch12$feedback.freq, ch12$PPI.person, do.exact = FALSE, do.asymp = TRUE)
#'
#' # Exercise 13.7 from "Applied Nonparametric Statistical Methods" (5th edition)
#' chisqtest.ANSM(ch13$medicine[ch13$location == "Rural"],
#'   ch13$response[ch13$location == "Rural"], seed = 1)
#'
#' @importFrom stats complete.cases chisq.test
#' @export
chisqtest.ANSM <-
  function(x, y = NULL, p = NULL, cont.corr = TRUE, max.exact.cases = 10,
           nsims.mc = 100000, seed = NULL, do.exact = TRUE, do.asymp = FALSE,
           do.mc = FALSE) {
    stopifnot((is.factor(x) && nlevels(x) > 1 | is.table(x)),
              (is.factor(y) | (is.table(x) && is.null(y))),
              (nlevels(y) > 1 | (is.table(x) && is.null(y))),
              (length(x) == length(y) | is.null(y)),
              ((is.table(x) && is.null(y) && is.numeric(p) &&
                  length(p) == length(x) && sum(p) == 1) | is.null(p)),
              is.logical(cont.corr),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.exact) == TRUE, is.logical(do.asymp) == TRUE,
              is.logical(do.mc) == TRUE)

    #labels
    varname1 <- deparse(substitute(x))
    if (!is.table(x)){
      varname2 <- deparse(substitute(y))
    }else{
      if (is.null(p)){
        varname2 <- NULL
      }else{
        varname2 <- deparse(substitute(p))
      }
    }

    #unused arguments
    H0 <- NULL
    alternative <- NULL
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
    if (is.table(x)){
      y <- NULL
      n <- sum(x)
    }else{
      complete.cases.id <- complete.cases(x, y)
      x <- x[complete.cases.id] #remove missing cases
      x <- droplevels(x)
      y <- y[complete.cases.id] #remove missing cases
      y <- droplevels(y)
      n <- length(x)
    }
    if (is.null(p)){
      probs <- rep(1/length(x), length(x))
    }else{
      probs <- p
    }
    suppressWarnings({
      chisq.test.out <- chisq.test(x, y, correct = cont.corr, p = probs)
    })
    stat <- chisq.test.out$statistic[[1]]

    #give Monte Carlo output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.mc <- TRUE
    }

    #exact p-value
    if(do.exact && n <= max.exact.cases){
      pval.exact.stat <- stat
      permutations <- perms(n)
      n.perms <- dim(permutations)[1]
      if (is.table(x)){
        x.vec <- NULL
        for (i in 1:length(x)){
          x.vec <- c(x.vec, rep(names(x)[i], x[i]))
        }
        x.vec <- factor(x.vec, levels = names(x))
      }else{
        x.vec <- x
      }
      pval.exact <- 0
      for (i in 1:n.perms){
        if (is.table(x)){
          suppressWarnings({
            chisq.tmp <- chisq.test(table(x.vec[permutations[i,]]), y,
                                    correct = cont.corr, p = probs)$statistic[[1]]
          })
        }else{
          suppressWarnings({
            chisq.tmp <- chisq.test(x[permutations[i,]], y,
                                    correct = cont.corr, p = probs)$statistic[[1]]
          })
        }
        if (chisq.tmp >= pval.exact.stat){
          pval.exact <- pval.exact + 1 / n.perms
        }
      }
    }

    #Monte Carlo p-value
    if (do.mc){
      pval.mc.stat <- stat
      if (!is.null(seed)){set.seed(seed)}
      suppressWarnings({
        pval.mc <- chisq.test(x, y, correct = cont.corr, p = probs,
                              simulate.p.value = TRUE, B = nsims.mc)$p.value
      })
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- stat
      pval.asymp <- chisq.test.out$p.value
    }

    #check if message needed
    if (!do.exact && !do.mc && !do.asymp) {
      test.note <- paste("Neither exact, asymptotic nor Monte Carlo test requested")
    }else if (do.exact && n > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact calculations\nrequired for ",
                          "exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ") so Monte ",
                          "Carlo p-value given")
    }

    #define hypotheses
    if (is.table(x)){
      if (is.null(p)){
        H0 <- paste0("H0: ", varname1, " follows a uniform distribution\n",
                     "H1: ", varname1, " does not follow a uniform distribution\n")
      }else{
        H0 <- paste0("H0: ", varname1, " follows the distribution defined by ",
                     varname2, "\n",
                     "H1: ", varname1, " does not follow the distribution ",
                     "defined by ", varname2, "\n")
      }
    }else{
      H0 <- paste0("H0: ", varname1, " and ", varname2, " are independent\n",
                   "H1: ", varname1, " and ", varname2, " are not independent\n")
    }

    #return
    result <- list(title = "Chi-squared test", varname1 = varname1,
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
