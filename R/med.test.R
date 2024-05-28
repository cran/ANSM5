#' Perform Median test
#'
#' @description
#' `med.test()` performs the Median test and is used in chapters 6 and 7 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector of same length as y
#' @param y Numeric vector, or factor of same length as x
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `1000`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @param do.CI Boolean indicating whether or not to perform confidence interval calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.7 from "Applied Nonparametric Statistical Methods" (5th edition)
#' med.test(ch6$males, ch6$females)
#'
#' # Example 7.5 from "Applied Nonparametric Statistical Methods" (5th edition)
#' med.test(ch7$time, ch7$surgeon, do.exact = FALSE, do.asymp = TRUE)
#'
#' @importFrom stats complete.cases median chisq.test
#' @export
med.test <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           CI.width = 0.95, max.exact.cases = 1000,
           do.asymp = FALSE, do.exact = TRUE, do.CI = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), (is.vector(y) && is.numeric(y)) |
              (is.factor(y) && length(x) == length(y) &&
                 length(x[complete.cases(x)]) == length(y[complete.cases(y)])),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              CI.width > 0, CI.width < 1, is.logical(do.asymp) == TRUE,
              is.logical(do.exact) == TRUE, is.logical(do.CI) == TRUE)
    alternative <- match.arg(alternative)

    #labels
    varname1 <- deparse(substitute(x))
    varname2 <- deparse(substitute(y))

    #default outputs
    cont.corr <- NULL
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
    y <- y[complete.cases(y)] #remove missing cases
    if (!is.null(H0)) {
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
    }
    if (!is.factor(y)){
      y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
      if (!is.null(H0)) {
        xy <- c(x - H0, y)
      }else{
        xy <- c(x, y)
      }
      med <- median(xy)
      x <- x[x != med] #ignoring x == med
      y <- y[y != med] #ignoring y == med
      nx <- length(x)
      ny <- length(y)
      n <- nx + ny
      x.gt <- sum(x > med)
      y.gt <- sum(y > med)
      propdiff <- x.gt / nx - y.gt / ny
      pbase <- factorial(nx) * factorial(ny) * factorial(x.gt + y.gt) *
        factorial(nx + ny - x.gt - y.gt) / factorial(nx + ny)
      tab <- table(c(rep("x", nx), rep("y", ny)), c(x, y) < med)
    }else{
      if (!is.null(H0)) {
        x2 <- x - H0
      }else{
        x2 <- x
      }
      med <- median(x2)
      x2 <- x2[x2 != med] #ignoring x2 == med
      n <- length(x2)
      g <- y
      table_g <- table(g)
      nlev <- nlevels(g)
      tab <- table(g, x2 < med)
      col1sum <- colSums(tab)[1]
      tab.p <- prod(factorial(colSums(tab))) * prod(factorial(rowSums(tab))) /
        (factorial(n) * prod(factorial(tab)))
    }

    #exact p-value
    if(!is.factor(y)){
      if (do.exact && n <= max.exact.cases){
        pval.exact <- 0
        for (x.gt.2 in 0:nx){
          y.gt.2 <- x.gt + y.gt - x.gt.2
          propdiffi <- x.gt.2 / nx - y.gt.2 / ny
          if (alternative == "less" && propdiffi <= propdiff){
            pval.exact <-
              pval.exact + pbase /
              (factorial(x.gt.2) * factorial(y.gt.2) * factorial(nx - x.gt.2) *
                 factorial(ny - y.gt.2))
          }else if (alternative == "greater" && propdiffi >= propdiff){
            pval.exact <-
              pval.exact + pbase /
              (factorial(x.gt.2) * factorial(y.gt.2) * factorial(nx - x.gt.2) *
                 factorial(ny - y.gt.2))
          }else if (alternative == "two.sided" &&
                    (abs(propdiffi) >= abs(propdiff) |
                     abs(abs(propdiffi) - abs(propdiff)) <
                     .Machine$double.eps ^ 0.5)){
            pval.exact <-
              pval.exact + pbase /
              (factorial(x.gt.2) * factorial(y.gt.2) * factorial(nx - x.gt.2) *
                 factorial(ny - y.gt.2))
          }
        }
      }
    }else{
      if (do.exact && n <= max.exact.cases){
        #evaluate all possible tables
        pval.exact <- 0
        ni <- rep(0, nlev)
        id <- nlev
        repeat{
          ni[id] <- ni[id] + 1
          if (ni[id] > table_g[nlev]){
            repeat{
              if (id == 2){break}
              ni[id - 1] <- ni[id - 1] + 1
              for (i in id:nlev){
                ni[i] <- 0
              }
              if (ni[id - 1] <= table_g[id - 1]){
                id <- nlev
                break
              }else{
                id <- id - 1
              }
            }
            if (id == 2){break}
          }
          ni[1] <- col1sum - sum(ni[2:nlev])
          if (ni[1] <= table_g[1] && ni[1] >= 0){
            tabi <- cbind(ni, table_g - ni)
            tabi.p <- prod(factorial(colSums(tabi))) *
              prod(factorial(rowSums(tabi))) /
              (factorial(n) * prod(factorial(tabi)))
            if (tabi.p <= tab.p){
              pval.exact <- pval.exact + tabi.p
            }
          }
        }
      }
    }

    #exact CI
    if (!is.factor(y) && do.CI && do.exact && n <= max.exact.cases){
      #prepare
      if (median(x) > median(y)){
        s1 <- x
        s2 <- y
      }else{
        s1 <- y
        s2 <- x
      }
      med <-median(c(s1, s2))
      s1.gt <- sum(s1 > med)
      s2.gt <- sum(s2 > med)
      ns1 <- s1.gt + sum(s1 < med) #ignoring s1 == med
      ns2 <- s2.gt + sum(s2 < med) #ignoring s2 == med
      #get increment
      increment <- min(abs(diff(sort(unique(c(s1, s2))))))
      ##get precision and round to avoid multiple dps created by differencing
      sorted <- sort(unique(c(s1, s2)))
      dp <- 0
      for (i in 1:length(sorted)){
        dp.pos <- unlist(gregexpr(".", as.character(sorted[i]), fixed = TRUE))
        if (dp.pos > 0 && nchar(as.character(sorted[i])) - dp.pos > dp){
          dp <- nchar(as.character(sorted[i])) - dp.pos
        }
      }
      increment <- round(increment, dp)
      #get range to try
      largest <- max(c(abs(s1), abs(s2)))
      #Lower limit
      for (i in seq(-largest, largest, increment)){
        s1a <- s1 + i
        median_tmp <- median(c(s1a, s2))
        a <- sum(s1a > median_tmp)
        b <- ns1 - a
        c <- s1.gt + s2.gt - a
        d <- ns2 - c
        mat <- matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE)
        pval.tmp <- fisher.test(mat)$p.value
        if (pval.tmp > (1 - CI.width)){break}
        CI.exact.lower <- i
        pval.lower <- pval.tmp / 2
      }
      #Upper limit
      for (i in seq(largest, CI.exact.lower + increment, -increment)){
        s1a <- s1 + i
        median_tmp <- median(c(s1a, s2))
        a <- sum(s1a > median_tmp)
        b <- ns1 - a
        c <- s1.gt + s2.gt - a
        d <- ns2 - c
        mat <- matrix(c(a, b, c, d), nrow = 2, ncol = 2, byrow = TRUE)
        pval.tmp <- fisher.test(mat)$p.value
        if (pval.tmp > (1 - CI.width)){
          CI.exact.upper <- i
          break
        }
        pval.upper <- pval.tmp / 2
      }
      #Actual CI width
      actualCIwidth.exact <- 1 - pval.lower - pval.upper
    }

    #asymptotic p-value
    if (do.asymp){
      asymp.test <- tryCatch(chisq.test(tab), warning=function(w)
        return(list(suppressWarnings(chisq.test(tab)), w)))
      if (length(asymp.test) ==2){
        pval.asymp.stat <- as.numeric(asymp.test[[1]]$statistic)
        pval.asymp <- asymp.test[[1]]$p.value
        pval.asymp.note <- paste0("NOTE: ", asymp.test[[2]]$message)
      }else{
        pval.asymp.stat <- as.numeric(asymp.test$statistic)
        pval.asymp <- asymp.test$p.value
        if (alternative != "two.sided"){
          pval.asymp <- pval.asymp / 2
        }
      }
    }
    #check if message needed
    if (!is.factor(y) && !do.asymp && !do.exact) {
      test.note <- paste("Neither exact nor asymptotic test/confidence ",
                         "interval requested")
    }else if (!do.asymp && !do.exact) {
        test.note <- paste("Neither exact nor asymptotic test requested")
    }else if (n > max.exact.cases) {
      affected <- NULL
      if (!is.factor(y) && do.exact && do.CI){
        affected <- "exact test and confidence interval"
      }else if (do.exact) {
        affected <- "exact test"
      }
      if (!is.null(affected)){
        test.note <- paste0("NOTE: Number of useful cases greater than ",
                            "current maximum allowed for exact\ncalculations ",
                            "required for ", affected, " (max.exact.cases = ",
                            sprintf("%1.0f", max.exact.cases), ")")
      }
    }
    #define hypotheses
    if (alternative == "two.sided"){
      H0 <- paste0("H0: samples are from populations with the same median\n",
                   "H1: samples are from populations with different medians\n")
    }else if (alternative == "less"){
      H0 <- paste0("H0: samples are from populations with the same median\n",
                   "H1: median of ", varname1, " is less than median of ",
                   varname2, "\n")
    }else if (alternative == "greater"){
      H0 <- paste0("H0: samples are from populations with the same median\n",
                   "H1: median of ", varname1, " is greater than median of ",
                   varname2, "\n")
    }
    #create title
    title <- "Median test"
    if (!is.factor(y) && do.exact && !is.null(pval.exact)){
      title <- "Median test (Fisher's Exact Test)"
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
