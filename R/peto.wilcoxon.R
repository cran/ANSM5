#' Perform Peto-Wilcoxon test
#'
#' @description
#' `peto.wilcoxon()` performs the Peto-Wilcoxon test and is used in chapter 9 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector of same length as y, x.c, y.c
#' @param y Numeric vector of same length as x, x.c, y.c
#' @param x.c Binary vector of same length as x, y, x.c
#' @param y.c Binary vector of same length as x, y, y.c
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `100000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 9.4 from "Applied Nonparametric Statistical Methods" (5th edition)
#' peto.wilcoxon(ch9$sampleI.survtime, ch9$sampleII.survtime,
#'   ch9$sampleI.censor, ch9$sampleII.censor, alternative = "less")
#'
#' @importFrom stats complete.cases
#' @importFrom utils combn
#' @export
peto.wilcoxon <-
  function(x, y, x.c, y.c, alternative=c("two.sided", "less", "greater"),
           max.exact.perms = 100000, nsims.mc = 10000, seed = NULL) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              is.vector(x.c), is.numeric(x.c), is.vector(y.c), is.numeric(y.c),
              all(x.c == 0 | x.c == 1), all(y.c == 0 | y.c == 1),
              length(x) == length(x.c), length(y) == length(y.c),
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed))
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #unused arguments
    cont.corr <- NULL
    CI.width <- NULL
    do.asymp <- FALSE
    do.exact <- TRUE
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
    x <- x[complete.cases(x, x.c)] #remove missing cases
    x.c <- x.c[complete.cases(x, x.c)] #remove missing cases
    y <- y[complete.cases(y, y.c)] #remove missing cases
    y.c <- y.c[complete.cases(y, y.c)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    xy <- c(x, y)
    xy.c <- c(x.c, y.c)
    n.x <- length(x)
    n.y <- length(y)
    n <- length(xy)
    n.perms <- choose(n, min(length(x), length(y)))
    #calculate statistics
    tab <- rbind(c(0, 0), table(xy, xy.c))
    tab <- cbind(tab, rep(NA, dim(tab)[1])) #column 3 of Table 9.4
    tab <- cbind(tab, rep(NA, dim(tab)[1])) #column 4 of Table 9.4
    tab <- cbind(tab, rep(NA, dim(tab)[1])) #column 5 of Table 9.4
    tab[1, 3] <- n
    tab[1, 4] <- 1
    tab[1, 5] <- n
    for (i in 2:dim(tab)[1]){
      tab[i, 3] <- tab[i - 1, 3] - tab[i, 1] - tab[i, 2]
      tab[i, 4] <- tab[i, 3] / n
      if(tab[i, 2] == 0 && tab[i - 1, 2] == 1){
        tab[i, 5] <- dim(tab)[1] - i
      }else{
        tab[i, 5] <- tab[i - 1, 5] - tab[i, 1]
      }
    }
    tab <- cbind(tab, rep(NA, dim(tab)[1])) #column 6 of Table 9.4
    tab[1, 6] <- 1
    prob <- 1
    denom <- n
    for (i in 2:dim(tab)[1]){
      if (tab[i, 1] > 0){
        tab[i, 6] <- prob * tab[i, 5] / denom
      }else{
        tab[i, 6] <- tab[i - 1, 6]
      }
      if (tab[i, 2] > 0){
        prob <- tab[i, 6]
        denom <- tab[i, 3]
      }
    }
    R <- array(NA, c(n, 4))
    R[, 1] <- tab[match(xy, row.names(tab)) - 1, 6]
    R[, 2] <- tab[match(xy, row.names(tab)), 6]
    R[, 2][xy.c == 1] <- 0
    R[, 3] <- R[, 1] + R[, 2]
    R[, 4] <- R[, 3] / 2
    W.vec <- n + 0.5 - n * R[, 4]
    W.x <- sum(W.vec[1:n.x])
    W.y <- sum(W.vec[(n.x + 1):n])
    if (alternative == "two.sided"){
      W <- max(W.x, W.y)
    }else if (alternative == "less"){
      W <- W.y
    }else{
      W <- W.x
    }

    #exact p-value
    OverflowState <- FALSE
    if (n.perms <= max.exact.perms){
      #try complete combinations
      try_result <- suppressWarnings(try(
        combins <- combn(n, min(length(x), length(y))), silent = TRUE)
      )
      if (any(class(try_result) == "try-error")){
        OverflowState <- TRUE
      }else{
        n.combins <- dim(combins)[2]
      }
    }
    if (n.perms > max.exact.perms | OverflowState){
      #use Monte Carlo
      if (!is.null(seed)){set.seed(seed)}
      n.combins <- nsims.mc
    }
    #evaluate all combinations
    tmp.pval <- 0
    for (i in 1:n.combins){
      #define data
      if (n.perms > max.exact.perms | OverflowState){
        tmp.combins <- sample(n, min(length(x), length(y)))
      }else{
        tmp.combins <- combins[,i]
      }
      tmp.x <- xy[tmp.combins]
      tmp.x.c <- xy.c[tmp.combins]
      tmp.y <- xy[-tmp.combins]
      tmp.y.c <- xy.c[-tmp.combins]
      tmp.xy <- c(tmp.x, tmp.y)
      tmp.xy.c <- c(tmp.x.c, tmp.y.c)
      #calculate statistics
      tmp.tab <- rbind(c(0, 0), table(tmp.xy, tmp.xy.c))
      tmp.tab <- cbind(tmp.tab, rep(NA, dim(tmp.tab)[1])) #column 3 of Table 9.4
      tmp.tab <- cbind(tmp.tab, rep(NA, dim(tmp.tab)[1])) #column 4 of Table 9.4
      tmp.tab <- cbind(tmp.tab, rep(NA, dim(tmp.tab)[1])) #column 5 of Table 9.4
      tmp.tab[1, 3] <- n
      tmp.tab[1, 4] <- 1
      tmp.tab[1, 5] <- n
      for (i in 2:dim(tmp.tab)[1]){
        tmp.tab[i, 3] <- tmp.tab[i - 1, 3] - tmp.tab[i, 1] - tmp.tab[i, 2]
        tmp.tab[i, 4] <- tmp.tab[i, 3] / n
        if(tmp.tab[i, 2] == 0 && tmp.tab[i - 1, 2] == 1){
          tmp.tab[i, 5] <- dim(tmp.tab)[1] - i
        }else{
          tmp.tab[i, 5] <- tmp.tab[i - 1, 5] - tmp.tab[i, 1]
        }
      }
      tmp.tab <- cbind(tmp.tab, rep(NA, dim(tmp.tab)[1])) #column 6 of Table 9.4
      tmp.tab[1, 6] <- 1
      tmp.prob <- 1
      tmp.denom <- n
      for (i in 2:dim(tmp.tab)[1]){
        if (tmp.tab[i, 1] > 0){
          tmp.tab[i, 6] <- tmp.prob * tmp.tab[i, 5] / tmp.denom
        }else{
          tmp.tab[i, 6] <- tmp.tab[i - 1, 6]
        }
        if (tmp.tab[i, 2] > 0){
          tmp.prob <- tmp.tab[i, 6]
          tmp.denom <- tmp.tab[i, 3]
        }
      }
      tmp.R <- array(NA, c(n, 4))
      tmp.R[, 1] <- tmp.tab[match(tmp.xy, row.names(tmp.tab)) - 1, 6]
      tmp.R[, 2] <- tmp.tab[match(tmp.xy, row.names(tmp.tab)), 6]
      tmp.R[, 2][tmp.xy.c == 1] <- 0
      tmp.R[, 3] <- tmp.R[, 1] + tmp.R[, 2]
      tmp.R[, 4] <- tmp.R[, 3] / 2
      tmp.W <- n + 0.5 - n * tmp.R[, 4]
      tmp.W.x <- sum(tmp.W[1:n.x])
      tmp.W.y <- sum(tmp.W[(n.x + 1):n])
      if (alternative == "two.sided"){
        tmp.W <- max(tmp.W.x, tmp.W.y)
      }else if (alternative == "less"){
        tmp.W <- tmp.W.y
      }else{
        tmp.W <- tmp.W.x
      }
      #check against actual test statistic
      if (tmp.W > W){
        tmp.pval <- tmp.pval + 1 / n.combins
      }
    }
    if (alternative == "two.sided"){
      tmp.pval <- min(1, tmp.pval * 2)
    }
    #output
    if (n.perms <= max.exact.perms && !OverflowState){
      pval.exact.stat <- W
      pval.exact <- tmp.pval
    }else{
      pval.mc.stat <- W
      pval.mc <- tmp.pval
    }

    #check if message needed
    if (n.perms > max.exact.perms) {
      test.note <- paste0("NOTE: Number of permutations required greater than ",
                          "current maximum allowed for exact calculations\n",
                          "required for exact test (max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ") so Monte ",
                          "Carlo p-value given")
    }else if (OverflowState){
      test.note <- paste0("NOTE: Insufficient memory for exact calculations",
                          "required for exact test\n(max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ") so Monte ",
                          "Carlo p-value given")
    }

    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: the survival times distributions are identical\n",
                   "H1: the survival times distribution of ", varname1,
                   " are different from those of ", varname2, "\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: the survival times distributions are identical\n",
                   "H1: the survival times distribution of ", varname1,
                   " are below those of ", varname2, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: the survival times distributions are identical\n",
                   "H1: the survival times distribution of ", varname1,
                   " are above those of ", varname2, "\n")
    }

    #return
    result <- list(title = "Peto-Wilcoxon test", varname1 = varname1,
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
