#' Calculate Cohen's kappa
#'
#' @description
#' `cohen.kappa()` calculates Cohen's kappa and is used in chapter 10 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y1 Factor of same length as y2, blocks and same levels as y2 and (if blocks not NULL) with 2 levels
#' @param y2 Factor of same length as y1, blocks and same levels as y1 and (if blocks not NULL) with 2 levels
#' @param blocks Factor of same length as y1, y2 or NULL (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `10`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `FALSE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMstat object with the results from applying the function
#' @examples
#' # Example 10.11 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cohen.kappa(ch10$dentistA, ch10$dentistB, do.asymp = TRUE, do.exact = FALSE,
#'   alternative = "greater")
#'
#' # Example 10.12 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cohen.kappa(ch10$questionnaire, ch10$demonstration, ch10$items)
#'
#' @importFrom stats complete.cases pnorm pchisq quantile
#' @export
cohen.kappa <-
  function(y1, y2, blocks = NULL,
           alternative = c("two.sided", "less", "greater"), CI.width = 0.95,
           max.exact.cases = 10, nsims.mc = 100000, seed = NULL,
           do.asymp = FALSE, do.exact = TRUE, do.CI = FALSE, do.mc = FALSE) {
    stopifnot(is.factor(y1), is.factor(y2), all(levels(y1) == levels(y2)),
              (nlevels(y1) == 2 | is.null(blocks)), length(y1) == length(y2),
              ((is.factor(blocks) && nlevels(blocks) > 1) | is.null(blocks)),
              (alternative == c("two.sided", "less", "greater") |
                 alternative == "two.sided" | is.null(blocks)),
              (length(y1) == length(blocks) | is.null(blocks)),
              is.numeric(CI.width), length(CI.width) == 1,
              CI.width > 0, CI.width < 1,
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE, is.logical(do.mc) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(y1))
    varname2 <- deparse(substitute(y2))
    if (is.null(blocks)){
      varname3 <- NULL
    }else{
      varname3 <- deparse(substitute(blocks))
    }

    #unused arguments
    cont.corr <- NULL
    #default outputs
    H0 <- NULL
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
    if (is.null(blocks)){
      complete.cases.ID <- complete.cases(y1, y2)
      blocks <- as.factor(rep("X", length(y1)))
      b <- 1
    }else{
      complete.cases.ID <- complete.cases(y1, y2, blocks)
      blocks <- blocks[complete.cases.ID] #remove missing cases
      b <- nlevels(blocks)
    }
    y1 <- y1[complete.cases.ID] #remove missing cases
    y2 <- y2[complete.cases.ID] #remove missing cases
    n <- length(y1)
    g <- nlevels(y1)

    #create statistics
    tab <- table(y1, y2, blocks)
    kappa <- rep(NA, b)
    sek <- rep(NA, b)
    for (i in 1:b){
      proptable <- prop.table(tab[, , i])
      po <- sum(diag(proptable))
      rSums <- rowSums(proptable)
      cSums <- colSums(proptable)
      pe <- sum(rSums * cSums)
      kappa[i] <- (po - pe) / (1 - pe)
      sek[i] <- sqrt(pe + pe ^ 2 - sum(rSums * cSums * (rSums + cSums))) /
        ((1 - pe) * sqrt(sum(tab[, , i])))
    }
    stat <- kappa
    if (b == 1){
      statlabel <- "Cohen's kappa"
    }else{
      statlabel <- paste0("Cohen's kappa for ", levels(blocks))
    }

    #give mc output if exact not possible
    if (do.exact && n > max.exact.cases){
      do.mc <- TRUE
    }

    #do asymptotic test when blocks present
    if ((do.exact | do.mc) && b > 1) {
      do.asymp <- TRUE
    }

    #exact p-value
    if(do.exact && b == 1 && n <= max.exact.cases){
      permutations <- perms(n)
      n.perms <- dim(permutations)[1]
      pval.exact <- 0
      for (i in 1:n.perms){
        tab.tmp <- table(y1[permutations[i,]], y2)
        proptable.tmp <- prop.table(tab.tmp)
        po.tmp <- sum(diag(proptable.tmp))
        rSums.tmp <- rowSums(proptable.tmp)
        cSums.tmp <- colSums(proptable.tmp)
        pe.tmp <- sum(rSums.tmp * cSums.tmp)
        kappa.tmp <- (po.tmp - pe.tmp) / (1 - pe.tmp)
        if (alternative == "two.sided"){
          if (abs(kappa.tmp) >= abs(kappa)){
            pval.exact <- pval.exact + 1 / n.perms
          }
        }else if (alternative == "less"){
          if (kappa.tmp <= kappa){
            pval.exact <- pval.exact + 1 / n.perms
          }
        }else if (alternative == "greater"){
          if (kappa.tmp >= kappa){
            pval.exact <- pval.exact + 1 / n.perms
          }
        }
      }
    }

    #Monte Carlo p-value
    if(do.mc && b == 1){
      if (!is.null(seed)){set.seed(seed)}
      pval.mc <- 0
      for (i in 1:nsims.mc){
        y1.sample <- y1[sample(n, n, replace = FALSE)]
        tab.tmp <- table(y1.sample, y2)
        proptable.tmp <- prop.table(tab.tmp)
        po.tmp <- sum(diag(proptable.tmp))
        rSums.tmp <- rowSums(proptable.tmp)
        cSums.tmp <- colSums(proptable.tmp)
        pe.tmp <- sum(rSums.tmp * cSums.tmp)
        kappa.tmp <- (po.tmp - pe.tmp) / (1 - pe.tmp)
        if (alternative == "two.sided"){
          if (abs(kappa.tmp) >= abs(kappa)){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }else if (alternative == "less"){
          if (kappa.tmp <= kappa){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }else if (alternative == "greater"){
          if (kappa.tmp >= kappa){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }
      }
    }

    #asymptotic p-value
    if(do.asymp){
      if (b == 1){
        pval.asymp.stat <- kappa / sek
        pval.asymp <- pnorm(pval.asymp.stat, lower.tail = FALSE)
      }else{
        kappa2 <- rep(NA, b)
        W <- rep(NA, b)
        for (i in 1:b){
          tab2 <- c(tab[1, 1, i], tab[1, 2, i] + tab[2, 1, i], tab[2, 2, i])
          N <- sum(tab2)
          pihat <- ((2 * tab2[1]) + tab2[2]) / (2 * N)
          kappa2[i] <- 1 - (tab2[2] / (2 * N * pihat * (1 - pihat)))
          V <- ((1 - kappa2[i]) / N) * ((1 - kappa2[i]) * (1 - 2 * kappa2[i]) +
                                          (kappa2[i] * (2 - kappa2[i])) / (2 * pihat * (1 - pihat)))
          W[i] <- 1 / V
        }
        kappabar <- sum(W * kappa2) / sum(W)
        pval.asymp.stat <- sum(W * (kappa2 - rep(kappabar, b)) ^ 2)
        pval.asymp <- pchisq(pval.asymp.stat, b - 1, lower.tail = FALSE)
      }
    }

    #confidence interval
    if (do.CI && b == 1){
      if (!is.null(seed)){set.seed(seed)}
      kappa.vec <- rep(NA, nsims.mc)
      for (i in 1:nsims.mc){
        mc.sample <- sample(seq(1, n), n, replace = TRUE)
        tab.tmp <- table(y1[mc.sample], y2[mc.sample])
        proptable.tmp <- prop.table(tab.tmp)
        po.tmp <- sum(diag(proptable.tmp))
        rSums.tmp <- rowSums(proptable.tmp)
        cSums.tmp <- colSums(proptable.tmp)
        pe.tmp <- sum(rSums.tmp * cSums.tmp)
        kappa.vec[i] <- (po.tmp - pe.tmp) / (1 - pe.tmp)
      }
      CI.mc.lower <- quantile(kappa.vec, (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      CI.mc.upper <- quantile(kappa.vec, 1 - (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      CI.mc.note <- "Confidence interval is basic bootstrap interval"
    }

    #check if message needed
    if ((do.exact | do.mc) && b > 1) {
      stat.note <- paste("Only asymptotic test available when blocks present")
    }else if (do.exact && n > max.exact.cases) {
      stat.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact calculations\nrequired for ",
                          "exact test (max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ") so Monte ",
                          "Carlo p-value given")
    }
    if (do.CI && b > 1){
      if (!is.null(stat.note)){
        stat.note <- paste0(stat.note, "\n")
      }
      stat.note <- paste0(stat.note, "NOTE: Confidence intervals not ",
                          "available when blocks present")
    }

    #create hypotheses
    if (do.exact | do.mc | do.asymp){
      if (b == 1){
        H0 <- paste0("H0: Cohen's kappa for ", varname1, " and ", varname2,
                     " is 0")
        if (alternative == "two.sided"){
          H0 <- paste0(H0, "\nH1: Cohen's kappa for ", varname1, " and ",
                       varname2, " is not 0")
        }else if(alternative == "greater"){
          H0 <- paste0(H0, "\nH1: Cohen's kappa for ", varname1, " and ",
                       varname2, " is greater than 0")
        }else{
          H0 <- paste0(H0, "\nH1: Cohen's kappa for ", varname1, " and ",
                       varname2, " is less than 0")
        }
        H0 <- paste0(H0, "\n")
      }else{
        H0 <- paste0("H0: Cohen's kappa for ", varname1, " and ", varname2,
                     "\nis the same for all groups in ", varname3)
        H0 <- paste0(H0, "\nH1: Cohen's kappa for ", varname1, " and ",
                     varname2, "\nis not the same for all groups in ", varname3)
        H0 <- paste0(H0, "\n")
      }
    }

    #return
    result <- list(title = "Cohen's kappa", varname1 = varname1,
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
