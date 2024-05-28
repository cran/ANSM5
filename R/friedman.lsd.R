#' Perform Least Significant Differences test after the Friedman test
#'
#' @description
#' `friedman.lsd()` performs the Least Significant Differences test after the Friedman test and is used in chapter 8 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Numeric vector of same length as groups, blocks
#' @param groups Factor of same length as y, blocks with levels such that length(y) == nlevels(groups) * nlevels(blocks)
#' @param blocks Factor of same length as y, groups with levels such that length(y) == nlevels(groups) * nlevels(blocks)
#' @param ids Vector of length 2 with elements both levels of groups
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 8.11 from "Applied Nonparametric Statistical Methods" (5th edition)
#' friedman.lsd(ch8$prey.preference, ch8$prey, ch8$larva, c("Cyclops", "Anopheles"))
#'
#' #  from "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @importFrom stats complete.cases pt
#' @export
friedman.lsd <-
  function(y, groups, blocks, ids) {
    stopifnot(is.vector(y), is.numeric(y), is.factor(groups), is.factor(blocks),
              length(y) == length(groups), length(groups) == length(blocks),
              length(y) == nlevels(groups) * nlevels(blocks),
              length(ids) == 2, all(ids %in% levels(groups)))

    #labels
    varname1 <- deparse(substitute(y))
    varname2 <- paste0(deparse(substitute(groups)), " (as groups) with ",
                       deparse(substitute(blocks)), " (as blocks)")

    #unused arguments
    alternative <- NULL
    cont.corr <- NULL
    CI.width <- NULL
    do.CI <- FALSE
    max.exact.cases <- NULL
    nsims.mc <- NULL
    do.asymp <- TRUE
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
    complete.cases.ID <- complete.cases(y, groups, blocks)
    y <- y[complete.cases.ID] #remove missing cases
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    groups <- groups[complete.cases.ID] #remove missing cases
    blocks <- blocks[complete.cases.ID] #remove missing cases
    b <- nlevels(blocks)
    g <- nlevels(groups)
    rank.tab <- simplify2array(by(y, blocks, rank, simplify = TRUE))
    Sr <- sum(rank.tab ** 2)
    St <- sum(rowSums(rank.tab) ** 2) / b
    C <- b * g * (g + 1) ** 2 / 4
    Tstat <- b * (g - 1) * (St - C) / (Sr - C)
    T1stat <- (b - 1) * (St - C) / (Sr - St)
    tot.rank1 <- sum(rank.tab[match(ids[1], levels(groups)),])
    tot.rank2 <- sum(rank.tab[match(ids[2], levels(groups)),])

    #check for ties
    tiesexist = !all(rank.tab == round(rank.tab,0)) # TRUE if ties exist in blocks

    #calculate lsd test statistic and asymptotic p-value
    pval.asymp.stat <- abs(tot.rank1 - tot.rank2) /
      sqrt(2 * b* (Sr - St) / ((b - 1) * (g - 1)))
    pval.asymp <- pt(pval.asymp.stat, (b - 1) * (g - 1), lower.tail = FALSE)

    #create hypotheses
    H0 <- paste0("H0: ", ids[1], " and ", ids[2], " are from the same population\n",
                 "H1: samples differ in location\n")

    #return
    result <- list(title = paste0("Least Significant Differences test after ",
                                  "Friedman test\n"),
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
