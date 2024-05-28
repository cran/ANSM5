#' Perform Cochran Q test
#'
#' @description
#' `cochran.q()` performs the Cochran Q test and is used in chapter 7 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Binary vector of same length as groups, blocks
#' @param groups Factor of same length as y, blocks with levels such that length(y) == nlevels(groups) * nlevels(blocks)
#' @param blocks Factor of same length as y, groups with levels such that length(y) == nlevels(groups) * nlevels(blocks)
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `100000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `100000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 7.8 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cochran.q(ch7$outcome, ch7$climb, ch7$member, do.exact = FALSE, do.asymp = TRUE)
#'
#' # Exercise 7.14 from "Applied Nonparametric Statistical Methods" (5th edition)
#' cochran.q(ch7$soc.media.use, ch7$participant, ch7$day, do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases pf pchisq
#' @export
cochran.q <-
  function(y, groups, blocks, max.exact.perms = 100000, nsims.mc = 100000,
           seed = NULL, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(y), is.numeric(y), all(y == 0 | y == 1),
              is.factor(groups), is.factor(blocks),
              length(y) == length(groups), length(groups) == length(blocks),
              length(y) == nlevels(groups) * nlevels(blocks),
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)

    #labels
    varname1 <- deparse(substitute(y))
    varname2 <- paste0(deparse(substitute(groups)), " (as groups) with ",
                       deparse(substitute(blocks)), " (as blocks)")

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
    y <- y[complete.cases(y)] #remove missing cases
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    b <- nlevels(blocks)
    g <- nlevels(groups)
    n.perms <- factorial(g) ** (b - 1)
    T0 <- simplify2array(by(y, groups, sum, simplify = TRUE))
    B0 <- simplify2array(by(y, blocks, sum, simplify = TRUE))
    N <- sum(y)
    Q <- (g * (g - 1) * sum(T0 ** 2) - (g - 1) * N ** 2) / (g * N - sum(B0 ** 2))

    #exact p-value
    if(do.exact && n.perms <= max.exact.perms){
      pval.exact.stat <- Q
      #create combinations for each row
      for (i in 1:g){
        if (i == 1){
          rowcomb <- matrix(1)
        }else{
          rowcomb.nrows <- nrow(rowcomb)
          tmp.mtx <- matrix(nrow = i * rowcomb.nrows, ncol = i)
          for(j in 1:i){
            tmp.mtx[(j - 1) * rowcomb.nrows + 1:rowcomb.nrows,] <-
              cbind(j, rowcomb + (rowcomb >= j))
          }
          rowcomb <- tmp.mtx
        }
      }
      rowcomb.nrows <- dim(rowcomb)[1]
      #loop around combinations
      pval.exact <- 0
      b_i <- b
      comb_i <- c(rep(1, b - 1), 0)
      repeat{
        comb_i[b_i] <- comb_i[b_i] + 1
        if (comb_i[b_i] > rowcomb.nrows){
          repeat{
            if (b_i == 2){break}
            comb_i[b_i - 1] <- comb_i[b_i - 1] + 1
            for (i in b_i:b){
              comb_i[i] <- 1
            }
            if (comb_i[b_i - 1] <= rowcomb.nrows){
              b_i <- b
              break
            }else{
              b_i <- b_i - 1
            }
          }
          if (b_i == 2){break}
        }
        tab_i <- NULL
        for (ib in 1:b){
          tab_i <- rbind(tab_i, y[blocks == levels(blocks)[ib]][rowcomb[comb_i[ib],]])
        }
        B_i <- apply(tab_i, 1, sum)
        T_i <- apply(tab_i, 2, sum)
        Q_i <- (g * (g - 1) * sum(T_i ** 2) - (g - 1) * N ** 2) /
          (g * N - sum(B_i ** 2))
        if (Q_i >= pval.exact.stat){
            pval.exact <- pval.exact + 1 / n.perms
        }
      }
    }

    #Monte Carlo p-value
    if (do.exact && n.perms > max.exact.perms){
      if (!is.null(seed)){set.seed(seed)}
      pval.mc.stat <- Q
      pval.mc <- 0
      for (i in 1:nsims.mc){
        tab_i <- NULL
        for (ib in 1:b){
          tab_i <- rbind(tab_i, sample(y[blocks == levels(blocks)[ib]], g))
        }
        B_i <- apply(tab_i, 1, sum)
        T_i <- apply(tab_i, 2, sum)
        Q_i <- (g * (g - 1) * sum(T_i ** 2) - (g - 1) * N ** 2) /
          (g * N - sum(B_i ** 2))
        if (Q_i >= pval.mc.stat){
          pval.mc <- pval.mc + 1 / nsims.mc
        }
      }
    }

    #asymptotic p-value
    if(do.asymp){
      pval.asymp.stat <- Q
      pval.asymp <- pchisq(Q, g - 1, lower.tail = FALSE)
    }

    #check if message needed
    if (!do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test requested")
    }else if (do.exact && n.perms > max.exact.perms) {
      test.note <- paste0("NOTE: Number of permutations required greater than ",
                          "current maximum allowed for exact calculations\n",
                          "required for exact test (max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ") so Monte ",
                          "Carlo p-value given")
    }

    #create hypotheses
    H0 <- paste0("H0: proportions of ", varname1, " are identical\n")
    H0 <- paste0(H0, "H1: proportions differ")
    H0 <- paste0(H0, "\n")

    #return
    result <- list(title = "Cochran Q test",
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
