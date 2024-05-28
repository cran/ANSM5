#' Create bootstrap confidence interval
#'
#' @description
#' `bs()` creates a bootstrap confidence interval and is used in chapter 14 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param x Numeric vector
#' @param y Numeric vector or NULL (defaults to `NULL`)
#' @param CI.width Confidence interval width (defaults to `0.95`)
#' @param nsims.bs Number of bootstrap samples to be taken (defaults to `10000`)
#' @param seed Random number seed to be used for Monte Carlo simulations (defaults to `NULL`)
#' @returns A list object object with the results from applying the function
#' @examples
#' # Example 14.5 from "Applied Nonparametric Statistical Methods" (5th edition)
#' bs(ch14$example14.2, nsims.bs = 2000, CI.width = 0.95, seed = 1)
#' bs(ch14$example14.2, nsims.bs = 2000, CI.width = 0.99, seed = 1)
#'
#' @importFrom stats median quantile

#' @export
bs <-

  function(x, y = NULL, CI.width = 0.95, nsims.bs = 10000, seed = NULL) {

    stopifnot(is.vector(x), is.numeric(x), is.vector(y) | is.null(y),

              is.numeric(y) | is.null(y), is.numeric(CI.width),

              length(CI.width) == 1, CI.width > 0, CI.width < 1,

              is.numeric(nsims.bs), length(nsims.bs) == 1,

              is.numeric(seed) | is.null(seed),

              length(seed) == 1 | is.null(seed))



    #Monte Carlo confidence interval

    if (is.null(y)){ #for one sample

      median.vec <- rep(NA, nsims.bs)

      n <- length(x)

      if (!is.null(seed)){set.seed(seed)}

      for (i in 1:nsims.bs){

        x.2 <- sample(x, n, replace = TRUE)

        median.vec[i] <- median(x.2)

      }

      bs.ci.lower <- quantile(median.vec, (1 - CI.width) / 2, na.rm = TRUE)

      bs.ci.upper <- quantile(median.vec, 1 - (1 - CI.width) / 2, na.rm = TRUE)

      bs.vec <- median.vec

    }else{ #for two samples

      mediandiff.vec <- rep(NA, nsims.bs)

      nx <- length(x)

      ny <- length(y)

      if (!is.null(seed)){set.seed(seed)}

      for (i in 1:nsims.bs){

        x.2 <- sample(x, nx, replace = TRUE)

        y.2 <- sample(y, ny, replace = TRUE)

        mediandiff.vec[i] <- median(x.2) - median(y.2)

      }

      bs.ci.lower <- quantile(mediandiff.vec, (1 - CI.width) / 2, na.rm = TRUE)

      bs.ci.upper <- quantile(mediandiff.vec, 1 - (1 - CI.width) / 2, na.rm = TRUE)

      bs.vec <- mediandiff.vec

    }



    #return

    return(list("vector" = bs.vec, "nsamples" = nsims.bs, "seed" = seed,
           "CI" = c(as.numeric(bs.ci.lower), as.numeric(bs.ci.upper)),
           "CI.width" = CI.width))
  }
