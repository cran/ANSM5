#' Perform Zelen test
#'
#' @description
#' `zelen()` performs the Zelen test and is used in chapter 13 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Binary factor of same length as y, z
#' @param y Binary factor of same length as x, z
#' @param z Factor of same length as x, y
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `1000000`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Section 13.2.5 from "Applied Nonparametric Statistical Methods" (5th edition)
#' zelen(ch13$drug, ch13$side.effects, ch13$age.group)
#'
#' # Example 13.3 from "Applied Nonparametric Statistical Methods" (5th edition)
#' zelen(ch13$machine, ch13$output.status, ch13$material.source)
#'
#' @importFrom stats complete.cases dhyper
#' @export
zelen <-
  function(x, y, z, max.exact.perms = 1000000, do.exact = TRUE) {
    stopifnot(is.factor(x), is.factor(y), is.factor(z),
              nlevels(x) == 2, nlevels(y) == 2, nlevels(z) > 1,
              length(x) == length(y), length(y) == length(z),
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.logical(do.exact) == TRUE)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))
    varname3 <- deparse(substitute(z))

    #unused arguments
    H0 <- NULL
    cont.corr <- NULL
    alternative <- NULL
    do.mc <- FALSE
    do.asymp <- FALSE
    do.CI <- FALSE
    nsims.mc <- NULL
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
    CI.width <- NULL
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
    nlevels.z <- nlevels(z)
    tab.zall <- table(x, y)
    tab.z <- array(0, dim = c(nlevels.z, 2, 2))
    rtot.z <- array(rep(NA, 2), dim = c(nlevels.z, 2))
    ctot.z <- array(rep(NA, 2), dim = c(nlevels.z, 2))
    for (i in 1:nlevels.z){
      tab.z[i,,] <- table(x[z == levels(z)[i]], y[z == levels(z)[i]])
      rtot.z[i,] <- rowSums(tab.z[i,,])
      ctot.z[i,] <- colSums(tab.z[i,,])
    }

    #work out permutations
    if (do.exact){
      max.zperm.pos <- rep(sum(tab.zall) + 1, nlevels.z)
      zperm.row <- rep(NA, nlevels.z)
      zperm.col <- rep(NA, nlevels.z)
      for (k in 1:nlevels.z){
        for (i in 1:2){
          for (j in 1:2){
            if (tab.z[k, i, j] == min(tab.z[k,,])){
              if (min(rtot.z[k,i], ctot.z[k,j]) < max.zperm.pos[k]){
                max.zperm.pos[k] <- min(rtot.z[k,i], ctot.z[k,j])
                zperm.row[k] <- i
                zperm.col[k] <- j
              }
            }
          }
        }
      }
      n.perms <- prod(max.zperm.pos)
    }else{
      n.perms <- 0
    }

    #exact p-value and confidence interval
    if((do.exact | do.CI) && n.perms <= max.exact.perms){
      dhyper.z <- NULL
      for (k in 1:nlevels.z){
        dhyper.z <- c(dhyper.z, dhyper(x = tab.z[k, 1, 1], m = ctot.z[k, 1],
                                       n = ctot.z[k, 2], k = rtot.z[k, 1]))
      }
      tab.tmp <- array(0, dim = c(nlevels.z, 2, 2))
      probs <- NULL
      zperm.pos <- rep(0, nlevels.z)
      id <- nlevels.z
      repeat{
        for (k in 1:nlevels.z){
          tab.tmp[k, zperm.row[k], zperm.col[k]] <- zperm.pos[k]
          tab.tmp[k, zperm.row[k], 3 - zperm.col[k]] <-
            rtot.z[k, zperm.row[k]] - tab.tmp[k, zperm.row[k], zperm.col[k]]
          tab.tmp[k, 3 - zperm.row[k], zperm.col[k]] <-
            ctot.z[k, zperm.col[k]] - tab.tmp[k, zperm.row[k], zperm.col[k]]
          tab.tmp[k, 3 - zperm.row[k], 3 - zperm.col[k]] <-
            rtot.z[k, 3- zperm.row[k]] - tab.tmp[k, 3 - zperm.row[k], zperm.col[k]]
        }
        #calculate probabilities
        dhyper.tmp <- NULL
        for (k in 1:nlevels.z){
          dhyper.tmp <- c(dhyper.tmp,
                          dhyper(x = tab.tmp[k, 1, 1], m = ctot.z[k, 1],
                                 n = ctot.z[k, 2], k = rtot.z[k, 1]))
        }
        if (sum(tab.tmp[, 1, 1]) == tab.zall[1, 1]){
          probs <- c(probs, prod(dhyper.tmp))
        }
        #update
        for (i in nlevels.z:1){
          zperm.pos[i] <- zperm.pos[i] + 1
          if (zperm.pos[i] > max.zperm.pos[i]){
            zperm.pos[i] <- 0
          }else{
            break
          }
        }
        if (sum(zperm.pos) == 0){break}
      }
      #compute statistic and p-value
      if (do.exact){
        pval.exact.stat <- prod(dhyper.z) / sum(probs)
        probs <- probs / sum(probs)
        pval.exact <- sum(probs[probs <= pval.exact.stat])
      }
    }

    #check if message needed
    if (!do.exact){
      test.note <- paste("Exact test not requested")
    }else if (do.exact && n.perms > max.exact.perms) {
      test.note <- paste0("NOTE: Number of permutations required greater than ",
                          "current maximum allowed for exact calculations\n",
                          "required for exact test (max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ")")
    }

    #define hypotheses
    H0 <- paste0("H0: Odds ratio for ", varname1, " by ", varname2,
                 " is the same for all levels of ", varname3, "\n",
                 "H1: Odds ratio for ", varname1, " by ", varname2,
                 " is not the same for all levels of ", varname3, "\n")

    #return
    result <- list(title = "Zelen test", varname1 = varname1,
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
