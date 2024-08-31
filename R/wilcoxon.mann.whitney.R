#' Perform Wilcoxon-Mann-Whitney test
#'
#' @description
#' `wilcoxon.mann.whitney()` performs the Wilcoxon-Mann-Whitney test and is used in chapters 6, 8, 9 and 12 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector, or factor with same levels as y
#' @param y Numeric vector, or factor with same levels as x
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Examples 6.1 and 6.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' wilcoxon.mann.whitney(ch6$groupA, ch6$groupB)
#'
#' # Exercise 12.4 from "Applied Nonparametric Statistical Methods" (5th edition)
#' wilcoxon.mann.whitney(ch12$feedback.satisfaction[ch12$PPI.person.2 == "Representative"],
#'   ch12$feedback.satisfaction[ch12$PPI.person.2 == "Researcher"],
#'   do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases wilcox.test
#' @importFrom utils combn tail
#' @export
wilcoxon.mann.whitney <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           cont.corr = TRUE, CI.width = 0.95, max.exact.cases = 1000,
           nsims.mc = 100000, seed = NULL, do.asymp = FALSE, do.exact = TRUE,
           do.mc = FALSE, do.CI = TRUE) {
    stopifnot((is.vector(x) && is.numeric(x) && is.vector(y) && is.numeric(y)) |
                (is.factor(x) && is.factor(y) && all(levels(x) == levels(y)) &&
                   is.null(H0)),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              is.logical(cont.corr) == TRUE, CI.width > 0, CI.width < 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

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
    y <- y[complete.cases(y)] #remove missing cases
    n.x <- length(x)
    n.y <- length(y)
    n <- n.x + n.y
    if (is.factor(x) && is.factor(y)){
      x <- as.numeric(x)
      y <- as.numeric(y)
    }else{
      x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    }
    if (!is.null(H0)) {
      xy <- c(x - H0, y)
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
      xy <- round(c(x, y), -floor(log10(sqrt(.Machine$double.eps))))
    }
    s <- rank(xy, ties.method = "average")
    ranksumx <- sum(s[1 : length(x)])
    ranksumy <- sum(s) - ranksumx
    mannwhitneyux <- ranksumx - length(x) * (length(x) + 1) / 2
    mannwhitneyuy <- ranksumy - length(y) * (length(y) + 1) / 2
    stat <- paste0("\n", ranksumx, " (rank sum from ", varname1, "), ",
                   ranksumy, " (rank sum from ", varname2, ")\n", mannwhitneyux,
                   " (Mann-Whitney U from ", varname1, "), ", mannwhitneyuy,
                   " (Mann-Whitney U from ", varname2, ")")

    #give MC output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.mc <- TRUE
    }

    #check for ties
    tiesexist = !all(s == round(s,0)) # TRUE if ties exist

    #exact p-value & CI
    OverflowState <- FALSE
    if (do.exact && tiesexist){
      try_result <- suppressWarnings(try(
        combins <- combn(length(xy), min(length(x), length(y))), silent = TRUE)
        )
      if (any(class(try_result) == "try-error")){
        OverflowState <- TRUE
      }
      if (OverflowState){
        do.mc <- TRUE
      }
    }
    if (do.exact && ((!tiesexist && n <= max.exact.cases) |
                     (tiesexist && !OverflowState))){
      pval.exact.stat <- stat
      if (!tiesexist){
        wilcox.test.output <- wilcox.test(x, y, alternative = alternative,
                                          mu = H0, exact = TRUE,
                                          conf.int = do.CI,
                                          conf.level = CI.width)
        pval.exact <- wilcox.test.output$p.value
        if (do.CI){
          CI.exact.lower <- wilcox.test.output$conf.int[1]
          CI.exact.upper <- wilcox.test.output$conf.int[2]
        }
      }else{
        nfrom <- min(length(x), length(y))
        permfrom <- s * 2
        permfrom <- sort(permfrom)
        permsums <- rep(0,sum(tail(permfrom, nfrom)))
        for (i in 1:dim(combins)[2]){
          tmpsum <- sum(permfrom[combins[,i]])
          permsums[tmpsum] <- permsums[tmpsum] + 1
        }
        ranksum <- ifelse(length(x) < length(y), ranksumx * 2, ranksumy * 2)
        if (length(x) < length(y) && alternative == "less"){
          pval.exact.less <-
            sum(permsums[1:ranksum]) / sum(permsums)
        }else if (length(x) < length(y) && alternative == "greater"){
          pval.exact.greater <-
            sum(permsums[ranksum:length(permsums)]) / sum(permsums)
        }else if (length(x) >= length(y) && alternative == "less"){
          pval.exact.less <-
            sum(permsums[ranksum:length(permsums)]) / sum(permsums)
        }else if (length(x) >= length(y) && alternative == "greater"){
          pval.exact.greater <-
            sum(permsums[1:ranksum]) / sum(permsums)
        }else{
          pval.exact.less <-
            sum(permsums[1:ranksum]) / sum(permsums)
          pval.exact.greater <-
            sum(permsums[ranksum:length(permsums)]) / sum(permsums)
        }
        if (alternative=="two.sided"){
          pval.exact <- min(1, min(pval.exact.less, pval.exact.greater) * 2)
        }else if (alternative == "less"){
          pval.exact <- pval.exact.less
        }else if (alternative == "greater"){
          pval.exact <- pval.exact.greater
        }
        if (do.CI){
          CI.exact.lower <- NULL
          CI.exact.upper <- NULL
        }
      }
    }

    #asymptotic p-value and CI (with/without continuity correction)
    if (do.asymp){
      wilcox.test.output <- wilcox.test(x, y, alternative = alternative,
                                        mu = H0, exact = FALSE,
                                        correct = cont.corr, conf.int = do.CI,
                                        conf.level = CI.width)
      pval.asymp.stat <- stat
      pval.asymp <- wilcox.test.output$p.value
      if (do.CI){
        if (alternative != "two.sided"){
          wilcox.test.output <- wilcox.test(x, y, alternative = "two.sided",
                                            mu = H0, exact = FALSE,
                                            correct = cont.corr, conf.int = TRUE,
                                            conf.level = CI.width)
        }
        CI.asymp.lower <- wilcox.test.output$conf.int[1]
        CI.asymp.upper <- wilcox.test.output$conf.int[2]
      }
    }

    #Monte Carlo p-value
    if(do.mc){
      pval.mc.stat <- stat
      stat.mc <- wilcox.test(x, y, exact = FALSE, correct = cont.corr,
                             conf.int = FALSE)$statistic
      if (!is.null(seed)){set.seed(seed)}
      pval.lt <- 0
      pval.gt <- 0
      for (i in 1:nsims.mc){
        xy.tmp <- sample(n, n, replace = FALSE)
        x.tmp <- xy[xy.tmp[1:n.x]]
        y.tmp <- xy[xy.tmp[(n.x + 1):n]]
        stat.tmp <-
          wilcox.test(x.tmp, y.tmp, exact = FALSE, correct = cont.corr,
                      conf.int = FALSE)$statistic
        if (stat.tmp <= stat.mc){
          pval.lt <- pval.lt + 1 / nsims.mc
        }
        if (stat.tmp >= stat.mc){
          pval.gt <- pval.gt + 1 / nsims.mc
        }
      }
      if (alternative == "two.sided"){
        pval.mc <- min(1, min(pval.lt, pval.gt) * 2)
      }else if (alternative == "less"){
        pval.mc <- pval.lt
      }else if (alternative == "greater"){
        pval.mc <- pval.gt
      }
    }

    #Bootstrap CI
    if(do.CI && (do.mc | (do.exact && is.null(CI.exact.lower)))){
      bs.ci.res <- bs(x = x, y = y, CI.width = CI.width, nsims.bs = nsims.mc,
                      seed = seed)$CI
      CI.mc.lower <- bs.ci.res[1]
      CI.mc.upper <- bs.ci.res[2]
      CI.mc.note <- paste0("Confidence interval for difference (", varname1,
                           " minus ", varname2, ")\nis basic bootstrap interval for the median")
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test/confidence interval ",
                         "requested")
    }else if (n > max.exact.cases) {
      affected <- NULL
      if (do.exact && do.CI){
        affected <- "exact test and confidence interval"
      }else if (do.exact) {
        affected <- "exact test"
      }
      if (!is.null(affected)){
        test.note <- paste0("NOTE: Number of useful cases greater than current ",
                            "maximum allowed for exact\ncalculations required ",
                            "for ", affected, " (max.exact.cases = ",
                            sprintf("%1.0f", max.exact.cases), ")")
      }
    }
    if (tiesexist && (n > max.exact.cases | OverflowState)){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      test.note <- paste0(test.note, "NOTE: Ties exist in data and sample ",
                          "too large for exact\ncalculations required ",
                          "for exact test")
    }
    if (tiesexist && !OverflowState && do.exact && do.CI){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      if (do.CI){
        test.note <- paste0(test.note, "NOTE: Ties exist in data so exact ",
                            "confidence interval\nnot available")
      }else if (!do.CI){
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic test")
      }
    }
    if (tiesexist && do.asymp){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      if (do.CI){
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic\ntest and confidence interval")
      }else if (!do.CI){
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic test")
      }
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: samples are from the same population\n",
                   "H1: samples differ in location\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: samples are from the same population\n",
                   "H1: location of ", varname1, " is less than location of ",
                   varname2, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: samples are from the same population\n",
                   "H1: location of ", varname1, " is greater than location of ",
                   varname2, "\n")
    }

    #return
    result <- list(title = "Wilcoxon-Mann-Whitney test", varname1 = varname1,
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
