% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pearson.R
\name{pearson}
\alias{pearson}
\title{Calculate Pearson correlation}
\usage{
pearson(
  x,
  y,
  alternative = c("two.sided", "less", "greater"),
  max.exact.cases = 10,
  nsims.mc = 1e+05,
  seed = NULL,
  do.asymp = FALSE,
  do.exact = TRUE,
  do.mc = FALSE
)
}
\arguments{
\item{x}{Numeric vector of same length as y}

\item{y}{Numeric vector of same length as x}

\item{alternative}{Type of alternative hypothesis (defaults to \code{two.sided})}

\item{max.exact.cases}{Maximum number of cases allowed for exact calculations (defaults to \code{10})}

\item{nsims.mc}{Number of Monte Carlo simulations to be performed (defaults to \code{100000})}

\item{seed}{Random number seed to be used for Monte Carlo simulations (defaults to \code{NULL})}

\item{do.asymp}{Boolean indicating whether or not to perform asymptotic calculations (defaults to \code{FALSE})}

\item{do.exact}{Boolean indicating whether or not to perform exact calculations (defaults to \code{TRUE})}

\item{do.mc}{Boolean indicating whether or not to perform Monte Carlo calculations (defaults to \code{FALSE})}
}
\value{
An ANSMstat object with the results from applying the function
}
\description{
\code{pearson()} calculates the Pearson correlation and is used in chapters 10 and 11 of "Applied Nonparametric Statistical Methods" (5th edition)
}
\examples{
# Section 10.1.2 from "Applied Nonparametric Statistical Methods" (5th edition)
pearson(ch10$q1, ch10$q2, alternative = "greater", do.asymp = TRUE, do.exact = FALSE)

# Example 11.2 from "Applied Nonparametric Statistical Methods" (5th edition)
pearson(ch11$parentlimit, ch11$reportedtime - 1 * ch11$parentlimit, alternative = "two.sided")

}
