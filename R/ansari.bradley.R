#' Perform Ansari-Bradley test
#'
#' @description
#' `ansari.bradley()` performs the Ansari-Bradley test and is used in chapter 6 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector
#' @param H0 Null hypothesis value (defaults to `NULL`)
#' @param alternative Type of alternative hypothesis (defaults to `two.sided`)
#' @param max.exact.cases Maximum number of cases allowed for exact calculations (defaults to `25`)
#' @param do.asymp Boolean indicating whether or not to perform asymptotic calculations (defaults to `FALSE`)
#' @param do.exact Boolean indicating whether or not to perform exact calculations (defaults to `TRUE`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Example 6.12 from "Applied Nonparametric Statistical Methods" (5th edition)
#' ansari.bradley(ch6$typeA, ch6$typeB)
#'
#' # Exercise 6.16 from "Applied Nonparametric Statistical Methods" (5th edition)
#' ansari.bradley(ch6$travel, ch6$politics)
#'
#' @importFrom stats complete.cases median ansari.test
#' @importFrom utils combn
#' @export
ansari.bradley <-
  function(x, y, H0 = NULL, alternative=c("two.sided", "less", "greater"),
           max.exact.cases = 25, do.asymp = FALSE, do.exact = TRUE) {
    stopifnot(is.vector(x), is.numeric(x), is.vector(y), is.numeric(y),
              ((is.numeric(H0) && length(H0) == 1) | is.null(H0)),
              is.numeric(max.exact.cases), length(max.exact.cases) == 1,
              is.logical(do.asymp) == TRUE, is.logical(do.exact) == TRUE)
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
    y <- y[complete.cases(y)] #remove missing cases
    x <- round(x, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    y <- round(y, -floor(log10(sqrt(.Machine$double.eps)))) #handle floating point issues
    #equalise medians
    x <- x + (median(y) - median(x))
    n.x <- length(x)
    n.y <- length(y)
    if (!is.null(H0)) {
      xy <- c(x - H0, y)
      varname1 <- paste0(varname1, " - ", H0)
    }else{
      H0 <- 0
      xy <- c(x, y)
    }
    n.xy <- length(xy)
    #allocate ranks and calculate statistics
    xyranks <- rank(xy, ties.method = "average")
    absd <- abs(xyranks - (n.xy + 1) / 2)
    absd.x <- sum(absd[1:n.x])
    absd.y <- sum(absd[(n.x + 1):n.xy])
    if (absd.x < absd.y){
      n.s <- n.x
      absd.s <- absd.x
    }else{
      n.s <- n.y
      absd.s <- absd.y
    }

    #give asymptotic output if exact not possible
    if (do.exact && n.xy > max.exact.cases){
      do.asymp <- TRUE
    }

    #exact p-value
    if (do.exact && n.xy <= max.exact.cases){
      if (alternative == "two.sided"){
        pval.exact.stat <- absd.s
        all.combn <- combn(n.xy, n.s)
        count <- 0
        for (i in 1:dim(all.combn)[2]){
          if (sum(absd[all.combn[,i]]) <= absd.s) {
            count <- count + 2
          }
        }
      }else if (alternative == "less"){
        pval.exact.stat <- absd.x
        all.combn <- combn(n.xy, n.x)
        count <- 0
        for (i in 1:dim(all.combn)[2]){
          if (sum(absd[all.combn[,i]]) <= absd.x) {
            count <- count + 1
          }
        }
      }else if (alternative == "greater"){
        pval.exact.stat <- absd.x
        all.combn <- combn(n.xy, n.x)
        count <- 0
        for (i in 1:dim(all.combn)[2]){
          if (sum(absd[all.combn[,i]]) >= absd.x) {
            count <- count + 1
          }
        }
      }
      pval.exact <- count / dim(all.combn)[2]
      if (pval.exact > 1){pval.exact <- 1}
    }

    #asymptotic p-value
    if (do.asymp){
      if (alternative == "two.sided"){
        pval.asymp.stat <- absd.s
      }else{
        pval.asymp.stat <- absd.x
      }
      res <- ansari.test(x = x, y = y, alternative=alternative, exact = FALSE,
                         conf.int = FALSE)
      pval.asymp <- res$p.value
      pval.asymp.note <-
        paste0("Asymptotic test statistic obtained using method shown in ANSM ",
               "book\nand p-value obtained using ansari.test function from ",
               "stats package")
    }

    #define hypotheses
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

    #check if message needed
    if (do.exact && n.xy > max.exact.cases) {
      test.note <- paste0("NOTE: Number of useful cases greater than current ",
                          "maximum allowed for exact\n calculations required ",
                          "(max.exact.cases = ",
                          sprintf("%1.0f", max.exact.cases), ")")
    }

    #return
    result <- list(title = "Ansari-Bradley test", varname1 = varname1,
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
