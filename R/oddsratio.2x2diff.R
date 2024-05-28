#' Perform test for difference in odds ratios
#'
#' @description
#' `oddsratio.2x2diff()` performs the test for difference in odds ratios and is used in chapter 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Binary factor of same length as y, z
#' @param y Binary factor of same length as x, z
#' @param z Binary factor of same length as x, y
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `1000000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 13.2 from "Applied Nonparametric Statistical Methods" (5th edition)
#' oddsratio.2x2diff(ch13$physical.activity, ch13$tv.viewing, ch13$gender,
#'   do.exact = FALSE, do.asymp = TRUE)
#' oddsratio.2x2diff(ch13$physical.activity, ch13$tv.viewing, ch13$gender,
#'   do.exact = FALSE, do.mc = TRUE, seed = 1, nsims = 10000)
#'
#' @importFrom stats complete.cases r2dtable quantile pnorm qnorm
#' @export
oddsratio.2x2diff <-
  function(x, y, z, alternative = c("two.sided", "less", "greater"),
           CI.width = 0.95, max.exact.perms = 1000000, nsims.mc = 100000,
           seed = NULL, do.exact = TRUE, do.asymp = FALSE, do.mc = FALSE,
           do.CI = TRUE) {
    stopifnot(is.factor(x), is.factor(y), is.factor(z),
              nlevels(x) == 2, nlevels(y) == 2, nlevels(z) == 2,
              length(x) == length(y), length(y) == length(z),
              length(CI.width) == 1, CI.width > 0, CI.width < 1,
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.exact) == TRUE, is.logical(do.asymp) == TRUE,
              is.logical(do.mc) == TRUE, is.logical(do.CI) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))
    varname3 <- deparse(substitute(z))

    #unused arguments
    H0 <- NULL
    cont.corr <- NULL
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
    complete.cases.id <- complete.cases(x, y, z)
    x <- x[complete.cases.id] #remove missing cases
    y <- y[complete.cases.id] #remove missing cases
    z <- z[complete.cases.id] #remove missing cases
    x <- droplevels(x)
    y <- droplevels(y)
    z <- droplevels(z)
    tab.z1 <- table(x[z == levels(z)[1]], y[z == levels(z)[1]])
    rtot.z1 <- rowSums(tab.z1)
    ctot.z1 <- colSums(tab.z1)
    if (min(tab.z1) == 0){
      phi.star.1 <- log(((tab.z1[1, 1] + 0.5) * (tab.z1[2, 2] + 0.5)) /
                          ((tab.z1[1, 2] + 0.5) * (tab.z1[2, 1] + 0.5)))
    }else{
      phi.star.1 <- log((tab.z1[1, 1] * tab.z1[2, 2]) /
                          (tab.z1[1, 2] * tab.z1[2, 1]))
    }
    var.phi.star.1 <- sum(1 / tab.z1)
    tab.z2 <- table(x[z == levels(z)[2]], y[z == levels(z)[2]])
    rtot.z2 <- rowSums(tab.z2)
    ctot.z2 <- colSums(tab.z2)
    if (min(tab.z1) == 0){
      phi.star.2 <- log(((tab.z2[1, 1] + 0.5) * (tab.z2[2, 2] + 0.5)) /
                          ((tab.z2[1, 2] + 0.5) * (tab.z2[2, 1] + 0.5)))
    }else{
      phi.star.2 <- log((tab.z2[1, 1] * tab.z2[2, 2]) /
                          (tab.z2[1, 2] * tab.z2[2, 1]))
    }
    var.phi.star.2 <- sum(1 / tab.z2)
    stat <- phi.star.1 - phi.star.2

    #give mc output if exact not possible
    if (do.exact){
      max.z1perm.pos <- max(tab.z1) + 1
      for (i in 1:2){
        for (j in 1:2){
          if (tab.z1[i, j] == min(tab.z1)){
            if (min(rtot.z1[i], ctot.z1[j]) < max.z1perm.pos){
              max.z1perm.pos <- min(rtot.z1[i], ctot.z1[j])
              z1perm.row <- i
              z1perm.col <- j
            }
          }
        }
      }
      max.z2perm.pos <- max(tab.z2) + 1
      for (i in 1:2){
        for (j in 1:2){
          if (tab.z2[i, j] == min(tab.z2)){
            if (min(rtot.z2[i], ctot.z2[j]) < max.z2perm.pos){
              max.z2perm.pos <- min(rtot.z2[i], ctot.z2[j])
              z2perm.row <- i
              z2perm.col <- j
            }
          }
        }
      }
      n.perms <- (max.z1perm.pos + 1) * (max.z2perm.pos + 1)
      if (n.perms > max.exact.perms){
        do.mc <- TRUE
      }
    }else{
      n.perms <- 0
    }

    #exact p-value
    if(do.exact && n.perms <= max.exact.perms){
      pval.exact.stat <- stat
      pval.exact <- 0
      stat.tmp <- rep(NA, n.perms)
      prob.tmp <- rep(NA, n.perms)
      counter <- 1
      tab.tmp1 <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
      tab.tmp2 <- matrix(c(0, 0, 0, 0), nrow = 2, ncol = 2)
      #loop for z1 table
      for (z1perm.pos in 0:max.z1perm.pos){
        #create z1 table
        tab.tmp1[z1perm.row, z1perm.col] <- z1perm.pos
        tab.tmp1[z1perm.row, 3 - z1perm.col] <-
          rtot.z1[z1perm.row] - tab.tmp1[z1perm.row, z1perm.col]
        tab.tmp1[3 - z1perm.row, z1perm.col] <-
          ctot.z1[z1perm.col] - tab.tmp1[z1perm.row, z1perm.col]
        tab.tmp1[3 - z1perm.row, 3 - z1perm.col] <-
          rtot.z1[3- z1perm.col] - tab.tmp1[3 - z1perm.row, z1perm.col]
        #calculate table probability
        dhyper.tmp1 <- dhyper(x = tab.tmp1[1, 1], m = ctot.z1[1],
                              n = ctot.z1[2], k = rtot.z1[1])
        #calculate phi
        if (min(tab.tmp1) == 0){
          tab.tmp1 <- tab.tmp1 + 0.5
        }
        phi.star.tmp1 <- log((tab.tmp1[1, 1] * tab.tmp1[2, 2]) /
                               (tab.tmp1[1, 2] * tab.tmp1[2, 1]))
        #loop for z2 table
        for (z2perm.pos in 0:max.z2perm.pos){
          #create z2 table
          tab.tmp2[z2perm.row, z2perm.col] <- z2perm.pos
          tab.tmp2[z2perm.row, 3 - z2perm.col] <-
            rtot.z2[z2perm.row] - tab.tmp2[z2perm.row, z2perm.col]
          tab.tmp2[3 - z2perm.row, z2perm.col] <-
            ctot.z2[z2perm.col] - tab.tmp2[z2perm.row, z2perm.col]
          tab.tmp2[3 - z2perm.row, 3 - z2perm.col] <-
            rtot.z2[3- z2perm.col] - tab.tmp2[3 - z2perm.row, z2perm.col]
          #calculate table probability
          dhyper.tmp2 <- dhyper(x = tab.tmp2[1, 1], m = ctot.z2[1],
                                n = ctot.z2[2], k = rtot.z2[1])
          #calculate phi
          if (min(tab.tmp2) == 0){
            tab.tmp2 <- tab.tmp2 + 0.5
          }
          phi.star.tmp2 <- log((tab.tmp2[1, 1] * tab.tmp2[2, 2]) /
                                 (tab.tmp2[1, 2] * tab.tmp2[2, 1]))
          #compare with observed statistic
          stat.tmp[counter] <- phi.star.tmp1 - phi.star.tmp2
          prob.tmp[counter] <- dhyper.tmp1 * dhyper.tmp2
          if (alternative == "less" && stat.tmp[counter] <= pval.exact.stat){
            pval.exact <- pval.exact + prob.tmp[counter]
          }else if (alternative == "greater" && stat.tmp[counter] >= pval.exact.stat){
            pval.exact <- pval.exact + prob.tmp[counter]
          }else if (alternative == "two.sided" &&
                    abs(stat.tmp[counter]) >= abs(pval.exact.stat)){
            pval.exact <- pval.exact + prob.tmp[counter]
          }
          counter <- counter + 1
        }
      }
      if (do.CI){
        prob.tmp2 <- cumsum(prob.tmp[order(stat.tmp)])
        stat.tmp2 <- sort(stat.tmp)
        i <- 1
        repeat{
          if (prob.tmp2[i] <= (1 - CI.width) / 2 &&
              prob.tmp2[i + 1] > (1 - CI.width) / 2) {
            CI.exact.lower <- stat.tmp2[i]
            break
          }
          i <- i + 1
        }
        prob.tmp2 <- cumsum(prob.tmp[order(stat.tmp, decreasing = TRUE)])
        stat.tmp2 <- sort(stat.tmp, decreasing = TRUE)
        i <- 1
        repeat{
          if (prob.tmp2[i] <= (1 - CI.width) / 2 &&
              prob.tmp2[i + 1] > (1 - CI.width) / 2) {
            CI.exact.upper <- stat.tmp2[i]
            break
          }
          i <- i + 1
        }
      }
    }

    #Monte Carlo p-value
    if (do.mc){
      pval.mc.stat <- stat
      if (!is.null(seed)){set.seed(seed)}
      pval.mc <- 0
      stat.tmp <- rep(NA, nsims.mc)
      for (i in 1:nsims.mc){
        #simulate tables
        tab.tmp1 <- r2dtable(1, rtot.z1, ctot.z1)[[1]]
        tab.tmp2 <- r2dtable(1, rtot.z2, ctot.z2)[[1]]
        #create z1 stat
        phi.star.tmp1 <- log((tab.tmp1[1, 1] * tab.tmp1[2, 2]) /
                               (tab.tmp1[1, 2] * tab.tmp1[2, 1]))
        #create z2 stat
        phi.star.tmp2 <- log((tab.tmp2[1, 1] * tab.tmp2[2, 2]) /
                               (tab.tmp2[1, 2] * tab.tmp2[2, 1]))
        #compare with observed statistic
        stat.tmp[i] <- phi.star.tmp1 - phi.star.tmp2
        if (alternative == "less" && stat.tmp[i] <= pval.mc.stat){
          pval.mc <- pval.mc + 1 / nsims.mc
        }else if (alternative == "greater" && stat.tmp[i] >= pval.mc.stat){
          pval.mc <- pval.mc + 1 / nsims.mc
        }else if (alternative == "two.sided" &&
                  abs(stat.tmp[i]) >= abs(pval.mc.stat)){
          pval.mc <- pval.mc + 1 / nsims.mc
        }
      }
      if (do.CI){
        CI.mc.lower <- quantile(stat.tmp, (1 - CI.width) / 2, na.rm = TRUE)[[1]]
        CI.mc.upper <- quantile(stat.tmp, 1 - (1 - CI.width) / 2, na.rm = TRUE)[[1]]
      }
    }

    #asymptotic p-value
    if (do.asymp){
      pval.asymp.stat <- (phi.star.1 - phi.star.2) /
        sqrt(var.phi.star.1 + var.phi.star.2)
      if (alternative == "less"){
        pval.asymp <- pnorm(pval.asymp.stat, lower.tail = TRUE)
      }else if (alternative == "greater"){
        pval.asymp <- pnorm(pval.asymp.stat, lower.tail = FALSE)
      }else if (alternative == "two.sided"){
        pval.asymp <- pnorm(abs(pval.asymp.stat), lower.tail = FALSE) * 2
      }
      #confidence intervals
      if (do.CI){
        CI.asymp.lower <- phi.star.1 - phi.star.2 -
          qnorm((1 - CI.width) / 2, lower.tail = FALSE) *
          sqrt(var.phi.star.1 + var.phi.star.2)
        CI.asymp.upper <- phi.star.1 - phi.star.2 +
          qnorm((1 - CI.width) / 2, lower.tail = FALSE) *
          sqrt(var.phi.star.1 + var.phi.star.2)
      }
    }

    #check if message needed
    if (!do.mc && !do.asymp) {
      test.note <- paste("Neither Asymptotic nor Monte Carlo test requested")
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: Odds ratio for ", varname1, " by ", varname2,
                   " is the same for level '", levels(z)[1], "' and level '",
                   levels(z)[2], "' of ", varname3, "\n",
                   "H1: Odds ratio for ", varname1, " by ", varname2,
                   " is not the same for level '", levels(z)[1], "' and level '",
                   levels(z)[2], "' of ", varname3, "\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: Odds ratio for ", varname1, " by ", varname2,
                   " is the same for level '", levels(z)[1], "' and level '",
                   levels(z)[2], "' of ", varname3, "\n",
                   "H1: Odds ratio for ", varname1, " by ", varname2,
                   " is smaller for level '", levels(z)[1], "' than level '",
                   levels(z)[2], "' of ", varname3, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: Odds ratio for ", varname1, " by ", varname2,
                   " is the same for level '", levels(z)[1], "' and level '",
                   levels(z)[2], "' of ", varname3, "\n",
                   "H1: Odds ratio for ", varname1, " by ", varname2,
                   " is larger for level '", levels(z)[1], "' than level '",
                   levels(z)[2], "' of ", varname3, "\n")
    }

    #return
    result <- list(title = "Test for difference in odds ratios",
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
