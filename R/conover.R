#' Perform Conover test using standard or squared ranks
#'
#' @description
#' `conover()` performs the Conover test using standard or squared ranks and is used in chapters 6 and 7 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector of same length as y
#' @param y Factor of same length as x
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param abs.ranks Boolean indicating whether absolute ranks to be used instead of squared ranks (defaults to `FALSE`)
#' @param max.exact.perms Maximum number of permutations allowed for exact calculations (defaults to `5000000`)
#' @param nsims.mc Number of Monte Carlo simulations to be performed (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.mc Boolean indicating whether or not to perform Monte Carlo calculations (defaults to `FALSE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.13 from "Applied Nonparametric Statistical Methods" (5th edition)
#' conover(ch6$typeA, ch6$typeB, do.exact = FALSE, do.asymp = TRUE)
#'
#' # Exercise 7.15 from "Applied Nonparametric Statistical Methods" (5th edition)
#' conover(ch7$braking.distance.2, ch7$initial.speed, do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases pnorm pchisq var
#' @importFrom utils combn
#' @export
conover <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           abs.ranks = FALSE, max.exact.perms = 5000000, nsims.mc = 10000,
           seed = NULL, do.asymp = FALSE, do.exact = TRUE, do.mc = FALSE) {
    stopifnot(is.vector(x), is.numeric(x), (is.vector(y) && is.numeric(y)) |
              (is.factor(y) && length(x) == length(y) &&
                 length(x[complete.cases(x)]) == length(y[complete.cases(y)])),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.logical(abs.ranks) == TRUE,
              is.numeric(max.exact.perms), length(max.exact.perms) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              length(seed) == 1 | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.mc) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #unused arguments
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
    x <- x[complete.cases(x)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- y[complete.cases(y)] #remove missing cases
    if (is.vector(y)){
      y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    }
    n.x <- length(x)
    n.y <- length(y)
    n.xy <- n.x + n.y
    if (!is.null(H0)) {
      x <- x - H0
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
    }
    if (!is.factor(y)){
      mean.x <- mean(x)
      mean.y <- mean(y)
      dev.x <- abs(x - mean.x)
      dev.y <- abs(y - mean.y)
      dev.xy <- c(dev.x, dev.y)
      if (abs.ranks){
        xyranks <- rank(dev.xy, ties.method = "average")
      }else{
        xyranks <- rank(dev.xy, ties.method = "average") ** 2
      }
      xyranks.x <- sum(xyranks[1:n.x])
      xyranks.y <- sum(xyranks[(n.x + 1):n.xy])
      if (xyranks.x < xyranks.y){
        n.s <- n.x
        xyranks.s <- xyranks.x
      }else{
        n.s <- n.y
        xyranks.s <- xyranks.y
      }
      n.perms <- choose(n.xy, n.s)
    }else{
      #ensure alternative hypothesis is "two.sided"
      if (alternative != "two.sided"){
        alternative <- "two.sided"
        test.note <- paste0("NOTE: As ", varname2, " is a factor, alternative ",
                            "hypothesis must be two-sided")
      }
      #first reorder factor with smallest groups first for exact calculation purposes
      y <- factor(y, levels = levels(y)[rank(table(y), ties.method = "random")])
      y.count <- table(y)
      n.perms <- 1
      tmp.n <- 0
      for (i in 1:(nlevels(y) - 1)){
        n.perms <- n.perms * choose(n.x - tmp.n, as.integer(y.count[i]))
        tmp.n <- tmp.n + as.integer(y.count[i])
      }
      mean.x <- simplify2array(by(x, y, mean, simplify = TRUE))
      dev.x <- abs(x - mean.x[as.integer(y)])
      rank.x <- rank(dev.x)
      tab.rank.x <- simplify2array(by(rank.x ** 2, y, sum, simplify = TRUE))
      Sk <- sum(tab.rank.x ** 2 / y.count)
      C <- sum(rank.x ** 2) ** 2 / n.x
      Sr <- sum(rank.x ** 4)
      T0 <- (n.x - 1) * (Sk - C) / (Sr - C)
    }

    #give mc output if exact not possible
    if (do.exact && n.perms > max.exact.perms){
      do.mc <- TRUE
    }

    #exact p-value
    if (do.exact && n.perms <= max.exact.perms){
      if (!is.factor(y)){
        if (alternative == "two.sided"){
          pval.exact.stat <- xyranks.s
          all.combn <- combn(n.xy, n.s)
          count <- 0
          for (i in 1:dim(all.combn)[2]){
            if (sum(xyranks[all.combn[,i]]) <= pval.exact.stat) {
              count <- count + 2
            }
          }
        }else if (alternative == "less"){
          pval.exact.stat <- xyranks.x
          all.combn <- combn(n.xy, n.x)
          count <- 0
          for (i in 1:dim(all.combn)[2]){
            if (sum(xyranks[all.combn[,i]]) <= pval.exact.stat) {
              count <- count + 1
            }
          }
        }else if (alternative == "greater"){
          pval.exact.stat <- xyranks.x
          all.combn <- combn(n.xy, n.x)
          count <- 0
          for (i in 1:dim(all.combn)[2]){
            if (sum(xyranks[all.combn[,i]]) >= pval.exact.stat) {
              count <- count + 1
            }
          }
        }
        pval.exact <- count / dim(all.combn)[2]
      }else{
        combins <- NULL
        for (ig in 1:(nlevels(y) - 1)){
          if (ig == 1){
            combins <- t(combn(n.x, y.count[ig]))
          }else{
            combins2 <- NULL
            for (i in 1:dim(combins)[1]){
              combins2 <- rbind(combins2,
                                cbind(matrix(rep(combins[i,],
                                                 choose(n.x - dim(combins)[2],
                                                        y.count[ig])),
                                             ncol = dim(combins)[2],
                                             byrow = TRUE),
                                      t(combn(setdiff(seq(1:n.x), combins[i,]),
                                              y.count[ig]))))
            }
            combins <- combins2
          }
        }
        n.combins <- dim(combins)[1]
        pval.exact.stat <- T0
        pval.exact <- 0
        for (i in 1:n.combins){
          combin_i <- c(combins[i, ], setdiff(seq(1:n.x), combins[i,]))
          x_i <- x[combin_i]
          mean.x_i <- simplify2array(by(x_i, y, mean, simplify = TRUE))
          dev.x_i <- abs(x_i - mean.x_i[as.integer(y)])
          rank.x_i <- rank(dev.x_i)
          tab.rank.x_i <- simplify2array(by(rank.x_i ** 2, y, sum, simplify = TRUE))
          Sk_i <- sum(tab.rank.x_i ** 2 / y.count)
          C_i <- sum(rank.x_i ** 2) ** 2 / n.x
          Sr_i <- sum(rank.x_i ** 4)
          T0_i <- (n.x - 1) * (Sk_i - C_i) / (Sr_i - C_i)
          if (T0_i >= pval.exact.stat){
            pval.exact <- pval.exact + 1
          }
        }
        pval.exact <- pval.exact / n.combins
      }
    }

    #Monte Carlo p-value
    if (do.mc){
      if (!is.null(seed)){set.seed(seed)}
      if (!is.factor(y)){
        if (alternative == "two.sided"){
          pval.mc.stat <- xyranks.s
        }else{
          pval.mc.stat <- xyranks.x
        }
        pval.mc <- 0
        for (i in 1:nsims.mc){
          xy.sim <- sample(c(x, y), n.xy, replace = FALSE)
          x.sim <- xy.sim[1:n.x]
          y.sim <- xy.sim[(n.x + 1):n.xy]
          mean.x.sim <- mean(x.sim)
          mean.y.sim <- mean(y.sim)
          dev.x.sim <- abs(x.sim - mean.x.sim)
          dev.y.sim <- abs(y.sim - mean.y.sim)
          dev.xy.sim <- c(dev.x.sim, dev.y.sim)
          xyranks.sim <- rank(dev.xy.sim, ties.method = "average") ** 2
          xyranks.x.sim <- sum(xyranks.sim[1:n.x])
          xyranks.y.sim <- sum(xyranks.sim[(n.x + 1):n.xy])
          xyranks.s.sim <- min(xyranks.x.sim, xyranks.y.sim)
          if (alternative == "two.sided"){
            if (xyranks.s.sim <= pval.mc.stat){
              pval.mc <- pval.mc + 1 / nsims.mc
            }
          }else if (alternative == "less"){
            if (xyranks.x.sim <= pval.mc.stat){
              pval.mc <- pval.mc + 0.5 / nsims.mc
            }
          }else if (alternative == "greater"){
            if (xyranks.x.sim >= pval.mc.stat){
              pval.mc <- pval.mc + 0.5 / nsims.mc
            }
          }
        }
      }else{
        pval.mc.stat <- T0
        pval.mc <- 0
        for (i in 1:nsims.mc){
          x.sim <- sample(x, n.x, replace = FALSE)
          mean.x_i <- simplify2array(by(x.sim, y, mean, simplify = TRUE))
          dev.x_i <- abs(x.sim - mean.x_i[as.integer(y)])
          rank.x_i <- rank(dev.x_i)
          tab.rank.x_i <- simplify2array(by(rank.x_i ** 2, y, sum, simplify = TRUE))
          Sk_i <- sum(tab.rank.x_i ** 2 / y.count)
          C_i <- sum(rank.x_i ** 2) ** 2 / n.x
          Sr_i <- sum(rank.x_i ** 4)
          T0_i <- (n.x - 1) * (Sk_i - C_i) / (Sr_i - C_i)
          if (T0_i >= pval.mc.stat){
            pval.mc <- pval.mc + 1 / nsims.mc
          }
        }
      }
    }

    #asymptotic p-value (https://stat.ethz.ch/pipermail/r-help/2004-March/047190.html)
    if (do.asymp){
      if (!is.factor(y)){
        if (alternative == "two.sided"){
          pval.asymp.stat <- xyranks.s
          test.mean <- n.s * mean(xyranks)
          test.var <- n.s * (1 - n.s / (n.xy - 1)) * (var(xyranks) *
                                                        (n.xy - 1) / n.xy)
          pval.asymp <- pnorm((pval.asymp.stat - test.mean) / sqrt(test.var),
                              lower.tail = TRUE) * 2
        }else{
          pval.asymp.stat <- xyranks.x
          test.mean <- n.x * mean(xyranks)
          test.var <- n.x * (1 - n.x / (n.xy - 1)) * (var(xyranks) *
                                                        (n.xy - 1) / n.xy)
          if (alternative == "greater"){
            pval.asymp <- pnorm((pval.asymp.stat - test.mean) / sqrt(test.var),
                                lower.tail = FALSE)
          }else if (alternative == "less"){
            pval.asymp <- pnorm((pval.asymp.stat - test.mean) / sqrt(test.var),
                                lower.tail = TRUE)
          }
        }
      }else{
        pval.asymp.stat <- T0
        pval.asymp <- pchisq(T0, nlevels(y) - 1,lower.tail = FALSE)
      }
    }

    #define hypotheses
    if (!is.factor(y)){
      if (alternative == "two.sided"){
        H0 <- paste0("H0: samples have the same variance\n",
                     "H1: samples have different variances\n")
      }else if (alternative == "less"){
        H0 <- paste0("H0: samples have the same variance\n",
                     "H1: variance of ", varname1, " is less than variance of ",
                     varname2, "\n")
      }else if (alternative == "greater"){
        H0 <- paste0("H0: samples have the same variance\n",
                     "H1: variance of ", varname1, " is greater than variance of ",
                     varname2, "\n")
      }
    }else{
      H0 <- paste0("H0: samples have the same variance\n",
                   "H1: samples have different variances\n")
    }

    #check if message needed
    if (do.exact && n.perms > max.exact.perms) {
      if (!is.null(test.note)){
        test.note <- paste0(test.note, "\n")
      }
      test.note <- paste0("NOTE: Number of permutations required greater than ",
                          "current maximum allowed\nfor exact calculations ",
                          "required for exact test (max.exact.perms = ",
                          sprintf("%1.0f", max.exact.perms), ")\nso Monte ",
                          "Carlo p-value given")
    }

    if (abs.ranks){
      title <- "Conover test using standard ranks"
    }else{
      title <- "Conover test using squared ranks"
    }

    #return
    result <- list(title = title, varname1 = varname1,
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
