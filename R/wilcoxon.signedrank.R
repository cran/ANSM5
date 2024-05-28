#' Perform Wilcoxon signed-rank test
#'
#' @description
#' `wilcoxon.signedrank()` performs the Wilcoxon signed-rank test and is used in chapters 3, 4 and 5 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param cont.corr Boolean indicating whether or not to use continuity correction (defaults to `TRUE`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 3.4 from "Applied Nonparametric Statistical Methods" (5th edition)
#' wilcoxon.signedrank(ch3$heartrates1, 70, "greater")
#'
#' # Exercise 5.12 from "Applied Nonparametric Statistical Methods" (5th edition)
#' wilcoxon.signedrank(ch5$kHz0.125 - ch5$kHz0.25, 0)
#'
#' @importFrom stats complete.cases pnorm qnorm
#' @export
wilcoxon.signedrank <-
  function(x, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           cont.corr = TRUE, CI.width = 0.95, max.exact.cases = 1000,
           do.asymp = FALSE, do.exact = TRUE, do.CI = TRUE) {
    stopifnot(is.vector(x), is.numeric(x),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(cont.corr) == TRUE, is.numeric(CI.width),
              length(CI.width) == 1, CI.width > 0, CI.width < 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))

    #default outputs
    varname2 <- NULL
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
    multiplier.test <- NULL
    multiplier.ci <- NULL

    #prepare for tests
    if (!is.null(H0)){
      s.test <- rank(abs(x - H0), ties.method = "average")
      s.test[x == H0] <- 0 #make ranks equal to zero for values equal to H0
      multiplier.test = 2 - all(s.test == round(s.test,0)) # multiplier of 2 if average ranks, otherwise multiplier of 1
      if (!is.null(H0)){
        ranksumplus <- sum(s.test[sign(x - H0) == 1]) * multiplier.test
        ranksumminus <- sum(s.test[sign(x - H0) == -1]) * multiplier.test
      }
    }

    #prepare for CIs
    if (do.CI){
      s.ci <- rank(abs(x), ties.method = "average")
      multiplier.ci = 2 - all(s.ci == round(s.ci,0)) # multiplier of 2 if average ranks, otherwise multiplier of 1
      Walsh.averages <- NULL
      for (i in 1:length(x)){
        for (j in i:length(x)){
          Walsh.averages <- c(Walsh.averages, (x[i] + x[j]) / 2)
        }
      }
    }

    #give asymptotic output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.asymp <- TRUE
    }

    #exact p-value
    if (!is.null(H0) && do.exact && n <= max.exact.cases){
      permfrom.test <- s.test * multiplier.test
      permfrom.test <- sort(permfrom.test)
      permsums.test <- rep(0, sum(permfrom.test) + 1)
      permsums.test[1] <- 1
      current.max <- 0
      for (i in 1:length(permfrom.test)){
        permsumsnow.test <- permsums.test
        to.add <- permfrom.test[i]
        j <- 1:(current.max + 1)
        k <- (j - 1) + to.add + 1
        gtzero <- permsumsnow.test[j] > 0
        permsums.test[k][gtzero] <- permsums.test[k][gtzero] + permsumsnow.test[j][gtzero]
        current.max <- current.max + to.add
      }
      permsums.test <- data.frame(0:(length(permsums.test) - 1), permsums.test)

      pval.exact.stat <- paste0(ranksumminus / multiplier.test,
                                " (sum of negative ranks), ",
                                ranksumplus / multiplier.test,
                                " (sum of positive ranks)")
      pval.exact.less <-
        sum(permsums.test[permsums.test[, 1] >= ranksumminus, 2]) / sum(permsums.test[, 2])
      pval.exact.greater <-
        sum(permsums.test[permsums.test[, 1] >= ranksumplus, 2]) / sum(permsums.test[, 2])
      if (alternative=="two.sided"){
        pval.exact <- min(1, min(pval.exact.less, pval.exact.greater) * 2)
      }else if (alternative == "less"){
        pval.exact <- pval.exact.less
      }else if (alternative == "greater"){
        pval.exact <- pval.exact.greater
      }
    }

    #exact CI
    if (do.exact && do.CI && n <= max.exact.cases){
      permfrom.ci <- s.ci * multiplier.ci
      permfrom.ci <- sort(permfrom.ci)
      permsums.ci <- rep(0, sum(permfrom.ci) + 1)
      permsums.ci[1] <- 1
      current.max <- 0
      for (i in 1:length(permfrom.ci)){
        permsumsnow.ci <- permsums.ci
        to.add <- permfrom.ci[i]
        j <- 1:(current.max + 1)
        k <- (j - 1) + to.add + 1
        gtzero <- permsumsnow.ci[j] > 0
        permsums.ci[k][gtzero] <- permsums.ci[k][gtzero] + permsumsnow.ci[j][gtzero]
        current.max <- current.max + to.add
      }
      permsums.ci <- data.frame(0:(length(permsums.ci) - 1), permsums.ci)

      i <- 0
      denomsum <- sum(permsums.ci[, 2])
      repeat{
        pval.tmp <-
          sum(permsums.ci[permsums.ci[, 1] <= i, 2]) / denomsum * 2
        if (pval.tmp > (1 - CI.width)) {break}
        i <- i + multiplier.ci
      }
      CI.exact.lower <- sort(Walsh.averages, decreasing = FALSE)[i / multiplier.ci]
      CI.exact.upper <- sort(Walsh.averages, decreasing = TRUE)[i / multiplier.ci]
      actualCIwidth.exact <- 1 - sum(permsums.ci[permsums.ci[, 1] <= i - multiplier.ci, 2]) /
        denomsum * 2
      if (is.null(CI.exact.lower) | is.null(CI.exact.upper) | is.null(actualCIwidth.exact)){
        actualCIwidth.exact <- NULL
        CI.exact.lower <- NULL
        CI.exact.upper <- NULL
      }
    }

    #asymptotic p-value (with/without continuity correction)
    if (!is.null(H0) && do.asymp){
      s.test.2 <- min(ranksumminus, ranksumplus) / multiplier.test
      #with ties
      if (multiplier.test == 2){
        pval.asymp.stat <- (s.test.2 - 0.5 * sum(abs(s.test))) / (0.5 * sqrt(sum(s.test ** 2)))
        pval.asymp.less <-
          pnorm((ranksumminus / 2 - 0.5 * sum(abs(s.test))) / (0.5 * sqrt(sum(s.test ** 2))),
                lower.tail = FALSE)
        pval.asymp.greater <-
          pnorm((ranksumplus / 2 - 0.5 * sum(abs(s.test))) / (0.5 * sqrt(sum(s.test ** 2))),
                lower.tail = FALSE)
      }else{
        #without ties
        pval.asymp.stat <- ((s.test.2 + 0.5 * (cont.corr == TRUE)) - n * (n + 1) / 4) /
          sqrt(n * (n + 1) * (2 * n + 1) / 24)
        pval.asymp.less <-
          pnorm(((ranksumminus + (0.5 - (ranksumminus > ranksumplus)) *
                    (cont.corr == TRUE)) - n * (n + 1) / 4) /
                  sqrt(n * (n + 1) * (2 * n + 1) / 24), lower.tail = FALSE)
        pval.asymp.greater <-
          pnorm(((ranksumplus + (0.5 - (ranksumminus < ranksumplus)) *
                    (cont.corr == TRUE)) - n * (n + 1) / 4) /
                  sqrt(n * (n + 1) * (2 * n + 1) / 24), lower.tail = FALSE)
      }
      if (alternative=="two.sided"){
        pval.asymp <- min(1, min(pval.asymp.less, pval.asymp.greater) * 2)
      }else if (alternative == "less"){
        pval.asymp <- pval.asymp.less
      }else if (alternative == "greater"){
        pval.asymp <- pval.asymp.greater
      }
      if (n < 20){
        pval.asymp.note <- paste("(WARNING: n is less than 20 so asymptotic",
                                 "test is not recommended)")
      }
    }

    #asymptotic CI
    if (do.asymp && do.CI){
      if (n <= max.exact.cases){
        if (multiplier.ci == 2){ #recalculate if multiplier of 2 used for exact test
          permfrom.ci <- 1:n
          permsums.ci <- rep(0,sum(permfrom.ci) + 1)
          permsums.ci[1] <- 1
          current.max <- 0
          for (i in 1:length(permfrom.ci)){
            permsumsnow.ci <- permsums.ci
            to.add <- permfrom.ci[i]
            j <- 1:(current.max + 1)
            k <- (j - 1) + to.add + 1
            gtzero <- permsumsnow.ci[j] > 0
            permsums.ci[k][gtzero] <- permsums.ci[k][gtzero] + permsumsnow.ci[j][gtzero]
            current.max <- current.max + to.add
          }
          permsums.ci <- data.frame(0:(length(permsums.ci) - 1), permsums.ci)
        }
        meanstat <- sum(permsums.ci[, 1] * permsums.ci[, 2]) / sum(permsums.ci[, 2])
        sdstat <- sqrt(sum(permsums.ci[, 1] ^ 2 * permsums.ci[, 2]) / sum(permsums.ci[, 2]) - meanstat ^ 2)
        S <- qnorm((1 - CI.width) / 2) * sdstat + meanstat
        CI.asymp.lower <- sort(Walsh.averages, decreasing = FALSE)[round(S, 0) + 1]
        CI.asymp.upper <- sort(Walsh.averages, decreasing = FALSE)[length(Walsh.averages) - round(S, 0)]
        if (n < 20){
          CI.asymp.note <- paste("(WARNING: n is less than 20 so asymptotic",
                                 "CI is not recommended)")
        }
      }else{
        CI.asymp.note <- paste0("NOTE: Number of useful cases greater than current ",
                                "maximum allowed for exact\ncalculations required ",
                                "for statistics for asymptotic confidence interval\n",
                                "(max.exact.cases = ",
                                sprintf("%1.0f", max.exact.cases), ")")
      }
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test/confidence interval ",
                         "requested")
    }else if (n > max.exact.cases) {
      affected <- NULL
      if (do.exact && do.CI && !is.null(H0)){
        affected <- "exact test and confidence interval"
      }else if (do.exact && do.CI && is.null(H0)) {
        affected <- "exact confidence interval"
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

    if (multiplier.test == 2 && do.asymp && !is.null(H0)){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      if (cont.corr){
        cont.corr <- FALSE #set to FALSE so passed appropriately to print
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic test and the\ncontinuity ",
                            "correction has not been used")
      }else{
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic test")
      }
    }
    if (multiplier.ci == 2 && do.asymp && do.CI){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      if (cont.corr){
        cont.corr <- FALSE #set to FALSE so passed appropriately to print
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic confidence interval\nand the ",
                            "continuity correction has not been used")
      }else{
        test.note <- paste0(test.note, "NOTE: Ties exist in data so mid-ranks ",
                            "used for asymptotic confidence interval")
      }
    }
    if (is.null(H0) && (do.exact | do.asymp)){
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      test.note <- paste0(test.note, "NOTE: No value for H0 given so no p-value ",
                          "can be calculated")
    }
    if (!is.null(H0) && (do.exact | do.asymp)){
      if (sum(x == H0) != 0){
        if (!is.null(test.note)){
          test.note <- paste0(test.note, "\n")
        }
        test.note <- paste0(test.note, "NOTE: At least one value is equal to ",
                            "the null hypothesis - method of calculating\nranks ",
                            "which includes values equal to the null hypothesis ",
                            "in the ranking before\nsetting their ",
                            "ranks to zero has been used")
      }
    }

    #return
    result <- list(title = "Wilcoxon signed-rank test", varname1 = varname1,
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
