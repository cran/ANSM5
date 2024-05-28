#' Calculate Spearman beta
#'
#' @description
#' `spearman.beta()` calculates the Spearman beta and is used in chapter 11 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param y Numeric vector of same length as x
#' @param x Numeric vector of same length as y
#' @param H0 Null hypothesis value (defaults to `NULL`)
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
#' # Example 11.3 from "Applied Nonparametric Statistical Methods" (5th edition)
#' spearman.beta(ch11$reportedtime, ch11$parentlimit, H0 = 1)
#' spearman.beta(ch11$reportedtime, ch11$parentlimit, H0 = 1, do.CI = TRUE)
#'
#' @importFrom stats lm complete.cases median approx qnorm
#' @export
spearman.beta <-
  function(y, x, H0 = NULL, alternative = c("two.sided", "less", "greater"),
           CI.width = 0.95, max.exact.cases = 10, nsims.mc = 100000,
           seed = NULL, do.asymp = FALSE, do.exact = TRUE, do.CI = FALSE,
           do.mc = FALSE) {
    stopifnot(is.numeric(y), is.numeric(x), length(y) == length(x),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(CI.width), length(CI.width) == 1,
              CI.width > 0, CI.width < 1,
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.numeric(nsims.mc), length(nsims.mc) == 1,
              is.numeric(seed) | is.null(seed),
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE,
              is.logical(do.CI) == TRUE, is.logical(do.mc) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- paste0(deparse(substitute(y)), " ~ ", deparse(substitute(x)))
    varname2 <- NULL
    varname3 <- NULL

    #unused arguments
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
    CI.sample.lower <- NULL
    CI.sample.upper <- NULL
    CI.sample.note <- NULL
    stat.note <- NULL

    #prepare
    complete.cases.id <- complete.cases(x, y)
    y <- y[complete.cases.id] #remove missing cases
    x <- x[complete.cases.id] #remove missing cases
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    x.c <- x - median(x)
    n <- length(y)

    #calculate estimate of beta
    outer.y <- outer(y, y, "-")
    outer.x.c <- outer(x.c, x.c, "-")
    bvals <- outer.y / outer.x.c
    bvals <- sort(bvals[upper.tri(bvals)])
    mid.b <- c(NA, bvals[1:(length(bvals) - 1)] + diff(bvals) / 2, NA)
    T <- sum(1:n)
    for (i in 2:(length(mid.b) - 1)){
      T[i] <- sum(x.c * rank(y - mid.b[i] * x.c))
    }
    T[length(mid.b)] <- -T[1]
    rs <- T / max(T)
    stat <- approx(T[!duplicated(T)], mid.b[!duplicated(T)], xout=0)$y
    statlabel <- "Spearman beta"

    if (!is.null(H0)){
      spearman.test <- spearman(x.c, y - H0 * x.c, alternative = alternative,
                                max.exact.cases = max.exact.cases,
                                nsims.mc = nsims.mc, seed = seed,
                                do.asymp = do.asymp, do.exact = do.exact,
                                do.mc = do.mc)
      pval <- spearman.test$pval
      pval.stat <- spearman.test$pval.stat
      pval.note <- spearman.test$pval.note
      pval.asymp <- spearman.test$pval.asymp
      pval.asymp.stat <- spearman.test$stat
      pval.asymp.note <- spearman.test$pval.asymp.note
      pval.exact <- spearman.test$pval.exact
      pval.exact.stat <- spearman.test$stat
      pval.exact.note <- spearman.test$pval.exact.note
      pval.mc <- spearman.test$pval.mc
      pval.mc.stat <- spearman.test$stat
      pval.mc.note <- spearman.test$pval.mc.note
      stat.note <- spearman.test$stat.note
    }

    #CI from exact/asymptotic distribution
    if (do.CI){
      if (n <= max.exact.cases){
        combins <- perms(n)
        n.perms <- dim(combins)[1]
        corrs <- rep(NA, n.perms)
        for (i in 1:n.perms){
          corrs[i] <- cor(combins[i,], 1:n, method = "spearman")
        }
        corrs.dist <- cumsum(table(corrs) / n.perms)
        corr.CI.limit <- as.numeric(names(corrs.dist[sum(corrs.dist <= (1 - CI.width) / 2)]))
        corr.CI.p <- corrs.dist[sum(corrs.dist <= (1 - CI.width) / 2)][[1]]
        CI.sample.upper <- approx(rs[!duplicated(rs)], c(NA, bvals)[!duplicated(rs)], xout = corr.CI.limit)$y
        CI.sample.lower <- approx(rs[!duplicated(rs)], c(bvals, NA)[!duplicated(rs)], xout = -corr.CI.limit)$y
        actualCIwidth.exact <- 1 - corr.CI.p * 2
      }else{
        z <- qnorm((1 - CI.width) / 2, lower.tail = FALSE)
        r.approx <- z / sqrt(n - 1)
        CI.sample.lower <- approx(rs[!duplicated(rs)], c(bvals, NA)[!duplicated(rs)], xout = r.approx)$y
        CI.sample.upper <- approx(rs[!duplicated(rs)], c(NA, bvals)[!duplicated(rs)], xout = -r.approx)$y
      }
    }

    #create hypotheses
    if (!is.null(H0)){
      H0val <- H0
      H0 <- paste0("H0: Spearman beta for ", varname1, " is ", H0val)
      if (alternative == "two.sided"){
        H0 <- paste0(H0, "\nH1: Spearman beta for ", varname1, " is not ",
                     H0val)
      }else if(alternative == "greater"){
        H0 <- paste0(H0, "\nH1: Spearman beta for ", varname1,
                     " is greater than ", H0val)
      }else{
        H0 <- paste0(H0, "\nH1: Spearman beta for ", varname1, " is less than ",
                     H0val)
      }
      H0 <- paste0(H0, "\n")
    }

    #return
    result <- list(title = "Spearman beta", varname1 = varname1,
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
