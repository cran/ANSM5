#' Calculate Noether approximation
#'
#' @description
#' `noether()` calculates the Noether approximation and is used in chapter 5 of "Applied Nonparametric Statistical Methods" (5th edition)
#'
#' @param p1 Probability (expressed as a number between 0 and 1)
#' @param alpha Level of significance (expressed as number between 0 and 1) (defaults to `0.05`)
#' @param power Power (expressed as number between 0 and 1) (defaults to `0.9`)
#' @returns An ANSMtest object with the results from applying the function
#' @examples
#' # Exercise 5.8 from "Applied Nonparametric Statistical Methods" (5th edition)
#' noether(p1 = 0.7534, alpha = 0.05, power = 0.9)
#'
#' # Exercise 5.16 from "Applied Nonparametric Statistical Methods" (5th edition)
#' noether(p1 = 0.8, alpha = 0.025, power = 0.9)
#'
#' @importFrom stats qnorm
#' @export
noether <-
  function(p1, alpha = 0.05, power = 0.90) {
    stopifnot(is.numeric(p1), length(p1) == 1, p1 > 0, p1 < 1,
              is.numeric(alpha), length(alpha) == 1, alpha > 0, alpha < 1,
              is.numeric(power), length(power) == 1, power > 0, power < 1)

    zalpha <- qnorm(alpha, lower.tail = FALSE)
    zbeta <- qnorm(1 - power, lower.tail = FALSE)
    n <- ((zalpha + zbeta) ^ 2) / (4 * (p1 - 0.5) ^ 2)

    return(n)
  }
