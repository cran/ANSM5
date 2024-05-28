#' Calculate Kendall's concordance
#'
#' @description
#' `kendall.concordance()` calculates Kendall's concordance and is used in chapter 10 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Numeric vector of same length as groups, blocks
#' @param groups Factor of same length as y, blocks with levels such that length(y) == nlevels(groups) * nlevels(blocks)
#' @param blocks Factor of same length as y, groups with levels such that length(y) == nlevels(groups) * nlevels(blocks)
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `100000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMstat object with the results from applying the function
#' @examples
#' # Exercise 10.11 from "Applied Nonparametric Statistical Methods" (5th edition)
#' kendall.concordance(ch10$marks, ch10$script, ch10$examiner, do.exact = FALSE, do.asymp = TRUE)
#' kendall.concordance(ch10$marks, ch10$examiner, ch10$script, do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases
#' @export
kendall.concordance <-
  function(y, groups, blocks, max.exact.perms = 100000, nsims.mc = 100000,
           seed = NULL, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(y), is.numeric(y), is.factor(groups), is.factor(blocks),
              length(y) == length(groups), length(groups) == length(blocks),
              length(y) == nlevels(groups) * nlevels(blocks),
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)

    #labels
    varname1 <- deparse(substitute(blocks))
    varname2 <- NULL
    varname3 <- NULL

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
    CI.sample.lower <- NULL
    CI.sample.upper <- NULL
    CI.sample.note <- NULL
    stat.note <- NULL

    #prepare
    complete.cases.ID <- complete.cases(y, groups, blocks)
    y <- y[complete.cases.ID] #remove missing cases
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    groups <- droplevels(groups[complete.cases.ID]) #remove missing cases
    blocks <- droplevels(blocks[complete.cases.ID]) #remove missing cases
    ##ensure full replicate
    full.blocks <- levels(blocks)[(rowSums(table(blocks, groups)) == nlevels(groups))]
    keep.blocks <- blocks %in% full.blocks
    y <- y[keep.blocks]
    groups <- droplevels(groups[keep.blocks])
    blocks <- droplevels(blocks[keep.blocks])
    ##calculate
    g <- nlevels(groups)
    b <- nlevels(blocks)
    n.perms <- factorial(g) ** (b - 1)
    friedman.stat <- friedman(y, groups, blocks, do.asymp = TRUE,
                              do.exact = FALSE)$pval.asymp.stat
    stat <- friedman.stat / (b * (g - 1))
    statlabel <- "Kendall's concordance"

    #exact p-value
    if(do.exact && n.perms <= max.exact.perms){
      friedman.out <- friedman(y, groups, blocks, do.exact = TRUE)
      pval.exact <- friedman.out$pval.exact
      pval.exact.note <- friedman.out$pval.exact.note
    }

    #Monte Carlo p-value
    if (do.exact && n.perms > max.exact.perms){
      friedman.out <- friedman(y, groups, blocks, seed = seed, do.exact = TRUE)
      pval.mc <- friedman.out$pval.mc
      pval.mc.note <- friedman.out$pval.mc.note
    }

    #asymptotic p-value
    if(do.asymp){
      friedman.out <- friedman(y, groups, blocks, do.asymp = TRUE,
                               do.exact = FALSE)
      pval.asymp <- friedman.out$pval.asymp
      pval.asymp.note <- friedman.out$pval.asymp.note
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      stat.note <- paste("Neither exact nor asymptotic test requested")
    }else if (do.exact && n.perms > max.exact.perms) {
      stat.note <- paste0("NOTE: Number of permutations required greater than ",
                          "current maximum allowed for exact calculations\n",
                          "required for exact test (max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ") so Monte ",
                          "Carlo p-value given")
    }

    #create hypotheses
    H0 <- paste0("H0: Kendall's concordance for ", varname1, " is 0")
    H0 <- paste0(H0, "\nH1: Kendall's concordance for ", varname1, " is not 0")
    H0 <- paste0(H0, "\n")

    #return
    result <- list(title = "Kendall's concordance", varname1 = varname1,
                   varname2 = varname2, varname3 = varname3, stat = stat,
                   statlabel = statlabel, H0 = H0,
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
                   CI.mc.note = CI.mc.note, CI.sample.lower = CI.sample.lower,
                   CI.sample.upper = CI.sample.upper, CI.sample.note = CI.sample.note,
                   stat.note = stat.note)
    class(result) <- "ANSMstat"
    return(result)
  }
